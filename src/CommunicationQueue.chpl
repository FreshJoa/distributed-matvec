module CommunicationQueue {

use CTypes;
use BlockDist;
use Time;

use FFI;
use ForeignTypes;
use ConcurrentAccessor;

config const communicationQueueBufferSize = 32000;
config const stagingBuffersBufferSize = 100;

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
                  basisStates : c_ptr(uint(64)), coeffs : c_ptr(coeffType), size : int) {
  assert(basisPtr != nil);
  assert(accessorPtr != nil);

  // count == 0 has to be handled separately because c_ptrTo(indices) fails
  // when the size of indices is 0.
  if size == 0 then return;

  // TODO: is the fact that we're allocating `indices` over and over again
  // okay performance-wise?
  // logDebug("ls_hs_state_index ...");
  var indices : [0 ..# size] int = noinit;
  ls_hs_state_index(basisPtr.deref().payload, size, basisStates, 1, c_ptrTo(indices), 1);
  ref accessor = accessorPtr.deref();
  foreach k in 0 ..# size {
    const i = indices[k];
    const c = coeffs[k];
    // Importantly, the user could have made a mistake and given us an
    // operator which does not respect the basis symmetries. Then we could
    // have that a |σ⟩ was generated that doesn't belong to our basis. In
    // this case, we should throw an error.
    if i >= 0 then accessor.localAdd(i, c);
              else halt("invalid index");
  }
}
// inline proc localProcess(const ref sigmas : [] uint(64),
//                          const ref coeffs : [] complex(128)) {
//   // logDebug("localProcess(" + sigmas:string + ", " + coeffs:string + ") ...");
//   assert(sigmas.size == coeffs.size);
//   localProcess(sigmas.size, c_const_ptrTo(sigmas), c_const_ptrTo(coeffs));
// }


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

  var _locks : [LocaleSpace] sync bool;
  // var _basisPtr : c_ptr(Basis);
  // var _accessorPtr : c_ptr(ConcurrentAccessor(eltType));
  // var numRemoteCalls : atomic int;
  var flushBufferTime : atomic real;
  var enqueueUnsafeTime : atomic real;
  var enqueueTime : atomic real;
  var localProcessTime : atomic real;

  proc init(type coeffType, ptrStore : [] (c_ptr(Basis), c_ptr(ConcurrentAccessor(coeffType)))) {
    this.coeffType = coeffType;
    this._remoteBuffers = [i in LocaleSpace] new RemoteBuffer(coeffType, _dom.size, i);
    this._bases = [i in LocaleSpace] ptrStore[i][0];
    this._accessors = [i in LocaleSpace] ptrStore[i][1];
    complete();
  }

  inline proc _lock(localeIdx : int) { _locks[localeIdx].writeEF(true); }
  inline proc _unlock(localeIdx : int) { _locks[localeIdx].readFE(); }

  proc _flushBuffer(localeIdx : int, release : bool = false) {
    var timer = new Timer();
    timer.start();

    ref size = _sizes[localeIdx];
    // logDebug("processOnRemote(" + localeIdx:string + ") ...");
    const mySize = size;
    if mySize == 0 {
      if release then _unlock(localeIdx);

      timer.stop();
      flushBufferTime.add(timer.elapsed());
      return;
    }

    ref remoteBuffer = _remoteBuffers[localeIdx];
    remoteBuffer.put(_basisStates[localeIdx], _coeffs[localeIdx], mySize);
    size = 0;

    if release then _unlock(localeIdx);
    const remoteBasis = _bases[localeIdx];
    const remoteAccessor = _accessors[localeIdx];
    const remoteBasisStates = remoteBuffer.basisStates;
    const remoteCoeffs = remoteBuffer.coeffs;
    on Locales[localeIdx] {
      // const basisStates : [0 ..# size] uint(64) = sigmas;
      // const coeffs : [0 ..# size] cs.eltType = cs;
      // copyComplete$.writeEF(true);
      // ref queue = globalAllQueues[here.id]; // :c_ptr(owned CommunicationQueue(eltType));
      // if queue == nil then
      //   halt("oops: queuePtr is null");
      // logDebug("Calling localProcess ...");
      // logDebug("  " + queuePtr.deref()._basisPtr:string);
      // queue!.numRemoteCalls.add(1);
      // queue!.localProcess(basisStates, coeffs);
      localProcess(remoteBasis, remoteAccessor, remoteBasisStates, remoteCoeffs, mySize);
    }
    // We need to wait for the remote copy to complete before we can reuse
    // `sigmas` and `cs`.
    // copyComplete$.readFF();
    // logDebug("end processOnRemote!");
    timer.stop();
    flushBufferTime.add(timer.elapsed());
  }

  inline proc _enqueueUnsafe(localeIdx : int,
                             count : int,
                             basisStates : c_ptr(uint(64)),
                             coeffs : c_ptr(coeffType),
                             release : bool) {
    var timer = new Timer();
    timer.start();

    ref offset = _sizes[localeIdx];
    assert(offset + count <= _dom.size);
    c_memcpy(c_ptrTo(_basisStates[localeIdx][offset]), basisStates, count:c_size_t * c_sizeof(uint(64)));
    c_memcpy(c_ptrTo(_coeffs[localeIdx][offset]), coeffs, count:c_size_t * c_sizeof(coeffType));
    offset += count;
    // So far, everything was done locally. Only `processOnRemote` involves
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
               in coeffs : c_ptr(coeffType)) {
    var timer = new Timer();
    timer.start();

    if localeIdx == here.id {
      var t2 = new Timer();
      t2.start();
      localProcess(_bases[localeIdx], _accessors[localeIdx], basisStates, coeffs, count);
      t2.stop();
      localProcessTime.add(t2.elapsed());

      timer.stop();
      enqueueTime.add(timer.elapsed());
      return;
    }

    while count > 0 {
      _lock(localeIdx);
      const remaining = min(_dom.size - _sizes[localeIdx], count);
      _enqueueUnsafe(localeIdx, remaining, basisStates, coeffs, release=false);
      _unlock(localeIdx);
      count -= remaining;
      basisStates += remaining;
      coeffs += remaining;
    }
    timer.stop();
    enqueueTime.add(timer.elapsed());
  }

  proc drain() {
    for localeIdx in LocaleSpace {
      // if _sizes[localeIdx] > 0 {
      _lock(localeIdx);
      _flushBuffer(localeIdx, release=true);
      // }
    }
    // Check
    foreach localeIdx in LocaleSpace do
      assert(_sizes[localeIdx] == 0);
  }
}

record RemoteBuffer {
  type coeffType;
  var size : int;
  var localeIdx : int;
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);

  proc postinit() {
    const rvf_size = size;
    on Locales[localeIdx] {
      basisStates = c_malloc(uint(64), rvf_size);
      coeffs = c_malloc(coeffType, rvf_size);
    }
  }

  proc put(localBasisStates : [] uint(64),
           localCoeffs : [] coeffType,
           size : int) {
    assert(size <= this.size);
    PUT(c_ptrTo(localBasisStates[0]), localeIdx, basisStates, size:c_size_t * c_sizeof(uint(64)));
    PUT(c_ptrTo(localCoeffs[0]), localeIdx, coeffs, size:c_size_t * c_sizeof(coeffType));
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

  proc flush(localeIdx : int) {
    _queue.deref().enqueue(localeIdx, _sizes[localeIdx],
                   c_const_ptrTo(_basisStates[localeIdx, 0]),
                   c_const_ptrTo(_coeffs[localeIdx, 0]));
    _sizes[localeIdx] = 0;
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
    // This function could potentially be vectorized, because
    // computing hashes should make it compute bound.
    for i in 0 ..# batchSize {
      const localeIdx = localeIdxOf(basisStatesPtr[i]);
      add(localeIdx, basisStatesPtr[i], coeffsPtr[i]:coeffType);
    }
  }
}

}
