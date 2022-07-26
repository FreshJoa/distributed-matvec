module CommunicationQueue {

use CTypes;
use BlockDist;
use Time;
use ChapelLocks;

use FFI;
use ForeignTypes;
use ConcurrentAccessor;
use StatesEnumeration;

config const communicationQueueBufferSize = 32000;
config const stagingBuffersBufferSize = 100;
config const stagingBuffersBlockSize = 128;

// inline proc getAddr(const ref p): c_ptr(p.type) {
//   // TODO can this use c_ptrTo?
//   return __primitive("_wide_get_addr", p): c_ptr(p.type);
// }

inline proc GET(addr, node, rAddr, size) {
  __primitive("chpl_comm_get", addr, node, rAddr, size);
}

inline proc PUT(addr, node, rAddr, size) {
  __primitive("chpl_comm_put", addr, node, rAddr, size);
}

proc localProcess(basisPtr : c_ptr(Basis), accessorPtr : c_ptr(ConcurrentAccessor(?coeffType)),
                  basisStates : c_ptr(uint(64)), coeffs : c_ptr(?t), size : int) {
  var timer = new Timer();
  timer.start();
  var allocateTimer = new Timer();
  var indexingTimer = new Timer();
  var accessTimer = new Timer();
  local {
    // count == 0 has to be handled separately because c_ptrTo(indices) fails
    // when the size of indices is 0.
    if size == 0 then return (0, 0, 0, 0);

    allocateTimer.start();
    var indices : [0 ..# size] int = noinit;
    allocateTimer.stop();

    indexingTimer.start();
    ls_hs_state_index(basisPtr.deref().payload, size, basisStates, 1, c_ptrTo(indices[0]), 1);
    indexingTimer.stop();

    accessTimer.start();
    ref accessor = accessorPtr.deref();
    foreach k in 0 ..# size {
      const i = indices[k];
      const c = coeffs[k]:coeffType;
      // Importantly, the user could have made a mistake and given us an
      // operator which does not respect the basis symmetries. Then we could
      // have that a |σ⟩ was generated that doesn't belong to our basis. In
      // this case, we should throw an error.
      if i >= 0 then accessor.localAdd(i, c);
                else halt("invalid index");
    }
    accessTimer.stop();
  }
  timer.stop();
  return (timer.elapsed(), allocateTimer.elapsed(),
          indexingTimer.elapsed(), accessTimer.elapsed());
}

// Handles the non-trivial communication of matrix elements between locales
class CommunicationQueue {
  type coeffType;
  const _dom : domain(1) = {0 ..# communicationQueueBufferSize};

  var _sizes : [LocaleSpace] int;
  var _basisStates : [LocaleSpace] [_dom] uint(64);
  var _coeffs : [LocaleSpace] [_dom] coeffType;
  var _remoteBuffers : [LocaleSpace] RemoteBuffer(coeffType);

  var _bases : [LocaleSpace] c_ptr(Basis);
  var _accessors : [LocaleSpace] c_ptr(ConcurrentAccessor(coeffType));

  // var _locks : [LocaleSpace] sync bool;
  var _locks : [LocaleSpace] chpl_LocalSpinlock;
  // var _basisPtr : c_ptr(Basis);
  // var _accessorPtr : c_ptr(ConcurrentAccessor(eltType));
  // var numRemoteCalls : atomic int;
  var flushBufferTime : atomic real;
  var flushBufferPutTime : atomic real;
  var flushBufferRemoteTime : atomic real;
  var enqueueUnsafeTime : atomic real;
  var enqueueTime : atomic real;
  var lockWaitingTimeHere : atomic real;
  var lockWaitingTimeRemote : atomic real;
  var localProcessTimeHere : atomic real;
  var localProcessTimeRemote : atomic real;

  proc init(type coeffType, ptrStore : [] (c_ptr(Basis), c_ptr(ConcurrentAccessor(coeffType)))) {
    this.coeffType = coeffType;
    this._remoteBuffers = [i in LocaleSpace] new RemoteBuffer(coeffType, _dom.size, i);
    this._bases = [i in LocaleSpace] ptrStore[i][0];
    this._accessors = [i in LocaleSpace] ptrStore[i][1];
    complete();
  }

  inline proc _lock(localeIdx : int) { _locks[localeIdx].lock(); }
  inline proc _unlock(localeIdx : int) { _locks[localeIdx].unlock(); }
  inline proc _tryLock(localeIdx : int) {
    ref l = _locks[localeIdx].l;
    return !(l.read() || l.testAndSet(memoryOrder.acquire));
  }
  // inline proc _lock(localeIdx : int) { _locks[localeIdx].writeEF(true); }
  // inline proc _unlock(localeIdx : int) { _locks[localeIdx].readFE(); }

  proc _flushBuffer(localeIdx : int, release : bool = false) {
    var timer = new Timer();
    var lockWaitingTimer = new Timer();
    var putTimer = new Timer();
    var remoteTimer = new Timer();
    timer.start();

    ref size = _sizes[localeIdx];
    const mySize = size;
    if mySize == 0 {
      if release then _unlock(localeIdx);

      timer.stop();
      flushBufferTime.add(timer.elapsed(), memoryOrder.relaxed);
      return;
    }

    ref remoteBuffer = _remoteBuffers[localeIdx];
    lockWaitingTimer.start();
    remoteBuffer.lock.lock();
    lockWaitingTimer.stop();
    putTimer.start();
    remoteBuffer.put(_basisStates[localeIdx], _coeffs[localeIdx], mySize);
    putTimer.stop();
    
    const remoteBasis = _bases[localeIdx];
    const remoteAccessor = _accessors[localeIdx];
    const remoteBasisStates = remoteBuffer.basisStates;
    const remoteCoeffs = remoteBuffer.coeffs;

    size = 0;
    if release then _unlock(localeIdx);

    remoteTimer.start();
    var localProcessTime : real;
    const localProcessTimePtr = c_ptrTo(localProcessTime);
    const mainLocaleIdx = here.id;
    on Locales[localeIdx] {
      var timer = new Timer();
      timer.start();
      localProcess(remoteBasis, remoteAccessor, remoteBasisStates, remoteCoeffs, mySize);
      timer.stop();
      const time = timer.elapsed();
      PUT(c_const_ptrTo(time), mainLocaleIdx, localProcessTimePtr, c_sizeof(real));
    }
    remoteTimer.stop();

    remoteBuffer.lock.unlock();
    timer.stop();

    flushBufferTime.add(timer.elapsed(), memoryOrder.relaxed);
    flushBufferPutTime.add(putTimer.elapsed(), memoryOrder.relaxed);
    flushBufferRemoteTime.add(remoteTimer.elapsed(), memoryOrder.relaxed);
    lockWaitingTimeRemote.add(lockWaitingTimer.elapsed(), memoryOrder.relaxed);
    localProcessTimeRemote.add(localProcessTime, memoryOrder.relaxed);
  }

  inline proc _enqueueUnsafe(localeIdx : int,
                             count : int,
                             basisStates : c_ptr(uint(64)),
                             coeffs : c_ptr(?t),
                             release : bool) {
    var timer = new Timer();
    timer.start();

    ref offset = _sizes[localeIdx];
    assert(offset + count <= _dom.size);
    c_memcpy(c_ptrTo(_basisStates[localeIdx][offset]), basisStates,
             count:c_size_t * c_sizeof(uint(64)));
    // c_memcpy(c_ptrTo(_coeffs[localeIdx][offset]), coeffs, count:c_size_t * c_sizeof(coeffType));
    const dst = c_ptrTo(_coeffs[localeIdx][offset]);
    foreach i in 0 ..# count do
      dst[i] = coeffs[i]:coeffType;
    offset += count;
    // So far, everything was done locally. Only `_flushBuffer` involves
    // communication.
    if offset == _dom.size then
      _flushBuffer(localeIdx, release);
    else
      if release then _unlock(localeIdx);

    timer.stop();
    enqueueUnsafeTime.add(timer.elapsed());
  }

  proc enqueue(localeIdx : int,
               in count : int,
               in basisStates : c_ptr(uint(64)),
               in coeffs : c_ptr(?t)) {
    var timer = new Timer();
    timer.start();

    if localeIdx == here.id {
      var t2 = new Timer();
      t2.start();
      localProcess(_bases[localeIdx], _accessors[localeIdx], basisStates, coeffs, count);
      t2.stop();
      localProcessTimeHere.add(t2.elapsed());

      timer.stop();
      enqueueTime.add(timer.elapsed());
      return;
    }

    var lockWaitingTimer = new Timer();
    while count > 0 {
      lockWaitingTimer.start();
      _lock(localeIdx);
      lockWaitingTimer.stop();
      const remaining = min(_dom.size - _sizes[localeIdx], count);
      _enqueueUnsafe(localeIdx, remaining, basisStates, coeffs, release=true);
      // _unlock(localeIdx);
      count -= remaining;
      basisStates += remaining;
      coeffs += remaining;
    }
    timer.stop();
    enqueueTime.add(timer.elapsed(), memoryOrder.relaxed);
    lockWaitingTimeHere.add(lockWaitingTimer.elapsed(), memoryOrder.relaxed);
  }

  proc tryEnqueue(localeIdx : int,
                  in count : int,
                  in basisStates : c_ptr(uint(64)),
                  in coeffs : c_ptr(?t)) {
    var timer = new Timer();
    timer.start();

    if localeIdx == here.id {
      var t2 = new Timer();
      t2.start();
      localProcess(_bases[localeIdx], _accessors[localeIdx], basisStates, coeffs, count);
      t2.stop();
      localProcessTimeHere.add(t2.elapsed());

      timer.stop();
      enqueueTime.add(timer.elapsed());
      return true;
    }

    // var lockWaitingTimer = new Timer();
    if _tryLock(localeIdx) {
      while count > 0 {
        // lockWaitingTimer.start();
        // _lock(localeIdx);
        // lockWaitingTimer.stop();
        const remaining = min(_dom.size - _sizes[localeIdx], count);
        _enqueueUnsafe(localeIdx, remaining, basisStates, coeffs, release=false);
        // _unlock(localeIdx);
        count -= remaining;
        basisStates += remaining;
        coeffs += remaining;
      }
      _unlock(localeIdx);
      timer.stop();
      enqueueTime.add(timer.elapsed(), memoryOrder.relaxed);
      return true;
    }
    else {
      timer.stop();
      enqueueTime.add(timer.elapsed(), memoryOrder.relaxed);
      return false;
    }
    // lockWaitingTimeHere.add(lockWaitingTimer.elapsed(), memoryOrder.relaxed);
  }

  proc drain() {
    var lockWaitingTimer = new Timer();
    for localeIdx in LocaleSpace {
      lockWaitingTimer.start();
      _lock(localeIdx);
      lockWaitingTimer.stop();
      _flushBuffer(localeIdx, release=true);
    }
    // Check
    foreach localeIdx in LocaleSpace do
      assert(_sizes[localeIdx] == 0);
    lockWaitingTimeHere.add(lockWaitingTimer.elapsed(), memoryOrder.relaxed);
  }
}

record RemoteBuffer {
  type coeffType;
  var size : int;
  var localeIdx : int;
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var lock : chpl_LocalSpinlock;

  // var localProcessTime : real;
  proc postinit() {
    const rvf_size = size;
    on Locales[localeIdx] {
      basisStates = c_malloc(uint(64), rvf_size);
      coeffs = c_malloc(coeffType, rvf_size);
    }
  }

  proc put(localBasisStates : c_ptr(uint(64)),
           localCoeffs : c_ptr(coeffType),
           size : int) {
    assert(size <= this.size);
    PUT(localBasisStates, localeIdx, basisStates, size:c_size_t * c_sizeof(uint(64)));
    PUT(localCoeffs, localeIdx, coeffs, size:c_size_t * c_sizeof(coeffType));
  }

  inline proc put(localBasisStates : [] uint(64),
                  localCoeffs : [] coeffType,
                  size : int) {
    put(c_ptrTo(localBasisStates[0]), c_ptrTo(localCoeffs[0]), size);
  }

  proc deinit() {
    if basisStates != nil {
      const rvf_basisStates = basisStates;
      const rvf_coeffs = coeffs;
      on Locales[localeIdx] {
        assert(rvf_basisStates != nil && rvf_coeffs != nil);
        c_free(rvf_basisStates);
        c_free(rvf_coeffs);
      }
    }
  }
}

record StagingBuffers {
  type coeffType;
  const _capacity : int = stagingBuffersBufferSize;
  var _sizes : [0 ..# numLocales] int;
  var _basisStates : [0 ..# numLocales, 0 ..# _capacity] uint(64);
  var _coeffs : [0 ..# numLocales, 0 ..# _capacity] coeffType;
  var _queue : c_ptr(owned CommunicationQueue(coeffType));
  var hashTime : real;
  var enqueueTime : real;
  var nonZeroCount : int;
  var totalCount : int;

  proc init(ref queue : owned CommunicationQueue(?t)) {
    this.coeffType = t;
    this._queue = c_ptrTo(queue);
  }

  proc init=(const ref other : StagingBuffers) {
    assert(other.locale == here); // we do not support assignment from remote
    this.coeffType = other.coeffType;
    this._capacity = other._capacity;
    this._queue = other._queue;
    complete();
    foreach i in LocaleSpace {
      const n = other._sizes[i];
      if n != 0 {
        this._sizes[i] = n;
        c_memcpy(c_ptrTo(this._basisStates[i, 0]), c_const_ptrTo(other._basisStates[i, 0]),
                 n:c_size_t * c_sizeof(uint(64)));
        c_memcpy(c_ptrTo(this._coeffs[i, 0]), c_const_ptrTo(other._coeffs[i, 0]),
                 n:c_size_t * c_sizeof(coeffType));
      }
    }
  }

  proc deinit() {
    flush();
  }

  // proc deinit() {
    // logDebug("StagingBuffers: spent ", hashTime, " in localeIdxOf (and related stuff) and ",
    //          enqueueTime, " in enqueue; nonZeroCount / totalCount = ", nonZeroCount:real / totalCount:real);
  // }

  proc flush(localeIdx : int) {
    var timer = new Timer();
    timer.start();
    _queue.deref().enqueue(localeIdx, _sizes[localeIdx],
                           c_const_ptrTo(_basisStates[localeIdx, 0]),
                           c_const_ptrTo(_coeffs[localeIdx, 0]));
    _sizes[localeIdx] = 0;
    timer.stop();
    enqueueTime += timer.elapsed();
  }
  proc flush() {
    foreach localeIdx in LocaleSpace do
      flush(localeIdx);
  }

  proc add(localeIdx : int, basisState : uint(64), coeff : coeffType) {
    ref n = _sizes[localeIdx];
    assert(n < _capacity);
    _basisStates[localeIdx, n] = basisState;
    _coeffs[localeIdx, n] = coeff;
    n += 1;
    if n == _capacity then
      flush(localeIdx);
  }
  proc add(batchSize : int,
           basisStatesPtr : c_ptr(uint(64)),
           coeffsPtr : c_ptr(?t)) {
    var timer = new Timer();

    if numLocales == 1 {
      // timer.start();
      _queue.deref().enqueue(here.id, batchSize, basisStatesPtr, coeffsPtr);
      // timer.stop();
      // enqueueTime += timer.elapsed();
      return;
    }

    param blockSize = 128;
    const numBlocks = batchSize / blockSize;
    const numRest = batchSize % blockSize;

    timer.start();
    var localeIdxs : c_array(int, blockSize);
    for i in 0 ..# numBlocks {
      foreach k in 0 ..# blockSize {
        localeIdxs[k] = localeIdxOf(basisStatesPtr[i * blockSize + k]);
      }
      for k in 0 ..# blockSize {
        const c = coeffsPtr[i * blockSize + k]:coeffType;
        // if c != 0 {
          timer.stop();
          add(localeIdxs[k], basisStatesPtr[i * blockSize + k], c);
          timer.start();
        // }
      }
    }
    for i in blockSize * numBlocks ..< batchSize {
      const c = coeffsPtr[i]:coeffType;
      // totalCount += 1;
      // if c != 0 {
        // nonZeroCount += 1;
      const localeIdx = localeIdxOf(basisStatesPtr[i]);
      timer.stop();
      add(localeIdx, basisStatesPtr[i], c);
      timer.start();
      // }
    }

    timer.stop();
    hashTime += timer.elapsed();
    // This function could potentially be vectorized, because
    // computing hashes should make it compute bound.
    // for i in 0 ..# batchSize {
    //   const localeIdx = localeIdxOf(basisStatesPtr[i]);
    //   add(localeIdx, basisStatesPtr[i], coeffsPtr[i]:coeffType);
    // }
  }
}

}
