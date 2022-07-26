module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;
use Time;
use DynamicIters;
use CommDiagnostics;
use ChapelLocks;

use FFI;
use ForeignTypes;
use ConcurrentAccessor;
use BatchedOperator;
use CommunicationQueue;

config const kVerboseComm : bool = false; 
config const kVerboseGetTiming : bool = false; 
config const kUseQueue : bool = true;

/* 
 */
private proc localDiagonalBatch(indices : range(int, BoundedRangeType.bounded, false),
                                matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                const ref representatives : [] uint(64)) {
  const batchSize = indices.size;
  // assert(workspace.size >= batchSize);
  // if workspace.size < batchSize then
  //   workspace.domain = {0 ..# batchSize};
  ls_internal_operator_apply_diag_x1(
    matrix.payload, batchSize, c_const_ptrTo(representatives[indices.low]),
    c_ptrTo(y[indices.low]), c_const_ptrTo(x[indices.low]));
  // ls_hs_operator_apply_diag_kernel(
  //   matrix.payload, batchSize,
  //   c_const_ptrTo(representatives[indices.low]), 1,
  //   c_ptrTo(workspace));
  // foreach i in indices {
  //   y[i] = x[i] * workspace[i - indices.low]:eltType;
  // }
}

config const matrixVectorDiagonalNumChunks : int = 10 * here.maxTaskPar;
config const matrixVectorOffDiagonalNumChunks : int = 32 * here.maxTaskPar;
config const matrixVectorMainLoopNumTasks : int = here.maxTaskPar;

private proc localDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                           const ref representatives : [] uint(64),
                           numChunks : int = min(matrixVectorDiagonalNumChunks,
                                                 representatives.size)) {
  const totalSize = representatives.size;
  const batchSize = (totalSize + numChunks - 1) / numChunks;
  var ranges : [0 ..# numChunks] range(int, BoundedRangeType.bounded, false) =
    chunks(0 ..# totalSize, numChunks);
  // var workspace : [0 ..# batchSize] complex(128) = noinit;
  forall r in ranges {
    localDiagonalBatch(r, matrix, x, y, representatives);
  }
}

var globalPtrStore : [LocaleSpace] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64))));

private proc localOffDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                              const ref representatives : [] uint(64),
                              numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                    representatives.size)) {
  var timer = new Timer();
  timer.start();
  var initTime : real;
  var computeOffDiagTime : atomic real;
  var stagingAddTime : atomic real;
  var stagingFlushTime : atomic real;
  var queueDrainTime : real;
  var initTimer = new Timer();
  initTimer.start();

  const chunkSize = (representatives.size + numChunks - 1) / numChunks;
  logDebug("Local dimension=", representatives.size, ", chunkSize=", chunkSize);
  var accessor = new ConcurrentAccessor(y);
  globalPtrStore[here.id] = (c_const_ptrTo(matrix.basis), c_ptrTo(accessor));
  allLocalesBarrier.barrier();

  var queue = new CommunicationQueue(eltType, globalPtrStore);

  const ranges : [0 ..# numChunks] range(int) =
    chunks(0 ..# representatives.size, numChunks);
  initTimer.stop();
  initTime += initTimer.elapsed();
  forall rangeIdx in dynamic(0 ..# numChunks, chunkSize=1,
                             numTasks=matrixVectorMainLoopNumTasks)
                     with (ref queue,
                           var batchedOperator = new BatchedOperator(matrix, chunkSize),
                           var staging = new StagingBuffers(queue)) {
    const r : range(int) = ranges[rangeIdx];
    // logDebug("Processing ", r, " ...");

    var timer = new Timer();
    timer.start();
    const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
        r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
    chpl_task_yield();
    var radixOffsets : c_array(int, 257);
    radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
    chpl_task_yield();
    // const n = offsetsPtr[r.size];
    timer.stop();
    computeOffDiagTime.add(timer.elapsed(), memoryOrder.relaxed);

    timer.clear();
    timer.start();
    var numberToDo = numLocales;
    var isComplete : [0 ..# numLocales] bool = false;
    while numberToDo > 0 {
      const prev = numberToDo;
      for localeIdx in 0 ..# numLocales {
        if isComplete[localeIdx] then continue;

        const k = radixOffsets[localeIdx];
        const n = radixOffsets[localeIdx + 1] - k;
        const isDone = if n != 0 then queue.tryEnqueue(localeIdx, n, basisStatesPtr + k, coeffsPtr + k)
                                 else true;
        if isDone {
          isComplete[localeIdx] = true;
          numberToDo -= 1;
        }
        // else
        //   chpl_task_yield();
      }
      if numberToDo == prev then chpl_task_yield();
    }
    // staging.add(n, basisStatesPtr, coeffsPtr);
    timer.stop();
    stagingAddTime.add(timer.elapsed(), memoryOrder.relaxed);
  }

  var queueTimer = new Timer();
  queueTimer.start();
  queue.drain();
  queueTimer.stop();
  queueDrainTime += queueTimer.elapsed();

  allLocalesBarrier.barrier();
  timer.stop();
  logDebug("localOffDiagonal took ", timer.elapsed(), "\n",
           "  ├─ ", initTime, " in initialization\n",
           "  ├─ ", computeOffDiagTime, " in computeOffDiag (total)\n",
           "  ├─ ", stagingAddTime, " in staging.add (total)\n",
           "  │   └─ ", queue.enqueueTime, " in queue.enqueue\n",
           "  │       ├─ ", queue.localProcessTimeHere, " in localProcess on here\n",
           "  │       ├─ ", queue.lockWaitingTimeHere, " waiting for queue.lock\n",
           "  │       └─ ", queue.enqueueUnsafeTime, " in queue._enqueueUnsafe\n",
           "  │           └─ ", queue.flushBufferTime, " in queue._flushBuffer\n",
           "  │               ├─ ", queue.lockWaitingTimeRemote, " waiting for remoteBuffer.lock\n",
           "  │               ├─ ", queue.flushBufferPutTime, " in remote PUTs\n",
           "  │               └─ ", queue.flushBufferRemoteTime, " in remote tasks\n",
           "  │                   └─ ", queue.localProcessTimeRemote, " in localProcess on remote\n",
           "  └─ ", queueDrainTime, " in queue.drain\n");
}

record PartitionInfo {
    var _countOrOffset : int;
    var nextOffset : int;

    inline proc count ref { return _countOrOffset; }
    inline proc offset ref { return _countOrOffset; }
};

proc partitionBy(in first : c_ptr(?eltType), last : c_ptr(eltType), predicate) {
    while true {
      if first == last then return last;
      if !predicate(first.deref()) then break;
      first += 1;
    }

    var it = first + 1;
    while it != last {
      if predicate(it.deref()) {
        first.deref() <=> it.deref();
        first += 1;
      }
      it += 1;
    }
    return first;
}

inline proc swapElements(a : int, b : int, arr : c_ptr(?t1)) {
  arr[a] <=> arr[b];
}
inline proc swapElements(a : int, b : int, arr1 : c_ptr(?t1), arr2 : c_ptr(?t2)) {
  swapElements(a, b, arr1);
  swapElements(a, b, arr2);
}

proc radixOneStep(numKeys : int, keys : c_ptr(uint(8)), offsets : c_array(int, 257), arrs...?numArrs)
{
  var partitions : c_array(PartitionInfo, 256);
  foreach i in 0 ..# numKeys {
    partitions[keys[i]:int].count += 1;
  }

  var remainingPartitions : c_array(uint(8), 256);
  var numPartitions : int;
  var total : int;
  for i in 0 ..# 256 {
    const count = partitions[i].count;
    if count > 0 {
      partitions[i].offset = total;
      total += count;
      remainingPartitions[numPartitions] = i:uint(8);
      numPartitions += 1;
    }
    partitions[i].nextOffset = total;
  }

  var lastRemaining = remainingPartitions:c_ptr(uint(8)) + numPartitions;
  var endPartition = remainingPartitions:c_ptr(uint(8)) + 1;
  while lastRemaining - endPartition > 0 {
    record Func {
      inline proc this(partitionIdx : uint(8)) {
        ref beginOffset = partitions[partitionIdx:int].offset;
        ref endOffset = partitions[partitionIdx:int].nextOffset;
        if beginOffset == endOffset then return false;

        for i in beginOffset .. endOffset - 1 {
          ref offset = partitions[keys[i]:int].offset;
          keys[i] <=> keys[offset];
          swapElements(i, offset, (...arrs));
          offset += 1;
        }
        return beginOffset != endOffset;
      }
    }
    lastRemaining = partitionBy(remainingPartitions:c_ptr(uint(8)), lastRemaining, new Func());
  }

  offsets[0] = 0;
  foreach i in 1 ..# 256 {
    offsets[i] = partitions[i - 1].nextOffset;
  }
}

record _LocalBuffer {
  type coeffType;
  var destLocaleIdx : int;
  var srcLocaleIdx : int;
  var capacity : int;
  var size : int;
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var isFull : c_ptr(chpl__processorAtomicType(bool));
  var isEmpty : chpl__processorAtomicType(bool);
  var isEOF : chpl__processorAtomicType(bool);

  proc postinit() {
    isEmpty.write(true);
    isEOF.write(false);
    assert(destLocaleIdx == here.id);
    // isFull.store(false);
    // const rvf_size = size;
    // on Locales[destLocaleIdx] {
      basisStates = c_malloc(uint(64), capacity);
      coeffs = c_malloc(coeffType, capacity);
      // isEmpty.deref().store(true);
    // }
  }

  proc deinit() {
    // assert(!isFull.read());
    // if basisStates != nil {
      // assert(coeffs != nil);
      // const rvf_basisStates = basisStates;
      // const rvf_coeffs = coeffs;
      // on Locales[destLocaleIdx] {
        c_free(basisStates);
        c_free(coeffs);
      // }
    // }
  }
}

record _RemoteBuffer {
  type coeffType;
  var destLocaleIdx : int;
  var srcLocaleIdx : int;
  var capacity : int;
  var size : c_ptr(int);
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var isFull : chpl__processorAtomicType(bool);
  var isEmpty : c_ptr(chpl__processorAtomicType(bool));
  var isEOF : c_ptr(chpl__processorAtomicType(bool));

  var _putTime : real;
  var _waitTime : real;
  var _submitTime : real;

  proc postinit() {
    isFull.write(false);
  }

  proc put(localBasisStates : c_ptr(uint(64)),
           localCoeffs : c_ptr(coeffType),
           count : int) {
    var timer = new Timer();
    timer.start();

    assert(here.id == srcLocaleIdx);
    assert(count <= capacity);
    assert(!isFull.read());
    PUT(localBasisStates, destLocaleIdx, basisStates, count:c_size_t * c_sizeof(uint(64)));
    PUT(localCoeffs, destLocaleIdx, coeffs, count:c_size_t * c_sizeof(coeffType));
    // PUT(c_const_ptrTo(count), destLocaleIdx, size, c_sizeof(int));

    timer.stop();
    _putTime += timer.elapsed();
  }

  inline proc put(localBasisStates : [] uint(64),
                  localCoeffs : [] coeffType,
                  size : int) {
    put(c_const_ptrTo(localBasisStates[0]), c_const_ptrTo(localCoeffs[0]), size);
  }

  proc submit(basisStatesPtr : c_ptr(uint(64)),
              coeffsPtr : c_ptr(coeffType),
              count : int,
              inTheMeantime) {
    var timer = new Timer();
    timer.start();

    var waitTimer = new Timer();
    waitTimer.start();
    // Wait for the buffer to become empty
    while isFull.read() {
      inTheMeantime();
      chpl_task_yield();
      // tryProcessLocal(taskIdx, newLocalBuffers, matrix, accessor);
    }
    waitTimer.stop();
    _waitTime += waitTimer.elapsed();

    put(basisStatesPtr, coeffsPtr, count);
    isFull.write(true);
    const atomicPtr = isEmpty;
    const sizePtr = size;
    on Locales[destLocaleIdx] {
      sizePtr.deref() = count;
      atomicStoreBool(atomicPtr, false);
    }

    timer.stop();
    _submitTime += timer.elapsed();
  }

  proc finish() {
    // startVerboseCommHere();
    const atomicPtr = isEOF;
    on Locales[destLocaleIdx] {
      atomicStoreBool(atomicPtr, true);
      // pragma "local fn" pragma "fast-on safe extern function"
      // extern proc atomic_store_bool(ref obj : chpl__processorAtomicType(bool), value : bool) : void;

      // atomic_store_bool(atomicPtr.deref(), true);
      // isEOF.deref().write(true);
    }
    // stopVerboseCommHere();
  }
}






var globalPtrStoreNoQueue : [LocaleSpace] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64))),
                                           c_ptr(_RemoteBuffer(complex(128))),
                                           c_ptr(_LocalBuffer(complex(128))));

config const remoteBufferSize = 150000;
config const numTasks = here.maxTaskPar;
config const kVerbose = false;
config const specialCase = false;

extern proc chpl_task_getId(): chpl_taskID_t;

inline proc atomicStoreBool(p : c_ptr(chpl__processorAtomicType(bool)), value : bool) {
  pragma "local fn" pragma "fast-on safe extern function"
  extern proc atomic_store_bool(ref _obj : chpl__processorAtomicType(bool),
                                _value : bool) : void;

  atomic_store_bool(p.deref(), value);
}

proc tryProcessLocal(taskIdx : int, srcLocaleIdx : int, newLocalBuffers, matrix, accessor) {
  ref localBuffer = newLocalBuffers[srcLocaleIdx, taskIdx];
  if !localBuffer.isEmpty.read() {
    // if kVerbose then
    //   logDebug("335: Calling localProcess for newLocalBuffers[", srcLocaleIdx, ", ", taskIdx, "]...");
    localProcess(c_const_ptrTo(matrix.basis), c_const_ptrTo(accessor),
                 localBuffer.basisStates,
                 localBuffer.coeffs,
                 localBuffer.size);
    localBuffer.isEmpty.write(true);
    const atomicPtr = localBuffer.isFull;
    on Locales[localBuffer.srcLocaleIdx] {
      atomicStoreBool(atomicPtr, false);
    }
    // if kVerbose then
    //   logDebug("335: done with localProcess for newLocalBuffers[", srcLocaleIdx, ", ", taskIdx, "]");
    return true;
  }
  return false;
}
proc tryProcessLocal(taskIdx : int, newLocalBuffers, matrix, accessor) {
  var hasDoneWork : bool = false;
  for srcLocaleIdx in 0 ..# numLocales {
    hasDoneWork = hasDoneWork ||
      tryProcessLocal(taskIdx, srcLocaleIdx, newLocalBuffers, matrix, accessor);
  }
  return hasDoneWork;
}

record LocalOffDiagonalTimers {
  var total : Timer;
  var initialization : Timer;
  var computeOffDiag : atomic real;
  var submit : real;
  var wait : real;
  var put : real;
}

proc _offDiagMakeLocalBuffers(numTasks : int, remoteBufferSize : int) {
  var localBuffers : [0 ..# numLocales, 0 ..# numTasks] _LocalBuffer(complex(128))
    = [(srcLocaleIdx, taskIdx) in {0 ..# numLocales, 0 ..# numTasks}]
        new _LocalBuffer(complex(128), here.id, srcLocaleIdx, remoteBufferSize);
  return localBuffers;
}

proc _offDiagMakeRemoteBuffers(numTasks : int, remoteBufferSize : int) {
  var remoteBuffers : [0 ..# numLocales, 0 ..# numTasks] _RemoteBuffer(complex(128))
    = [(destLocaleIdx, taskIdx) in {0 ..# numLocales, 0 ..# numTasks}]
        new _RemoteBuffer(complex(128), destLocaleIdx, here.id, remoteBufferSize);
  return remoteBuffers;
}

proc _offDiagInitLocalBuffers(ref localBuffers, const ref ptrStore) {
  const destLocaleIdx = localBuffers.locale.id;
  coforall loc in Locales do on loc {
    const srcLocaleIdx = loc.id;
    const (_basis, _accessor, _remoteBufferPtr, _localBufferPtr) = ptrStore[srcLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref remoteBuffer = (_remoteBufferPtr + destLocaleIdx * numTasks + taskIdx).deref();
      localBuffers[srcLocaleIdx, taskIdx].isFull = c_ptrTo(remoteBuffer.isFull);
    }
  }
}

proc _offDiagInitRemoteBuffers(ref remoteBuffers, const ref ptrStore) {
  const srcLocaleIdx = remoteBuffers.locale.id;
  coforall loc in Locales do on loc {
    const destLocaleIdx = loc.id;
    const (_basis, _accessor, _remoteBufferPtr, _localBufferPtr) = ptrStore[destLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref myLocalBuffer = (_localBufferPtr + srcLocaleIdx * numTasks + taskIdx).deref();
      ref remoteBuffer = remoteBuffers[destLocaleIdx, taskIdx];
      remoteBuffer.basisStates = myLocalBuffer.basisStates;
      remoteBuffer.coeffs = myLocalBuffer.coeffs;
      remoteBuffer.size = c_ptrTo(myLocalBuffer.size);
      remoteBuffer.isEmpty = c_ptrTo(myLocalBuffer.isEmpty);
      remoteBuffer.isEOF = c_ptrTo(myLocalBuffer.isEOF);
    }
  }
}

private proc localOffDiagonalNoQueue(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                     const ref representatives : [] uint(64),
                                     numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                           representatives.size)) {
  var timers = new LocalOffDiagonalTimers();
  timers.total.start();
  timers.initialization.start();

  // Other locales need access to our newLocalBuffers, newRemoteBuffers, etc.
  // so we create and store them in a common location (on locale 0) before a
  // barrier.
  var newLocalBuffers = _offDiagMakeLocalBuffers(numTasks, remoteBufferSize);
  var newRemoteBuffers = _offDiagMakeRemoteBuffers(numTasks, remoteBufferSize);
  var accessor = new ConcurrentAccessor(y);
  globalPtrStoreNoQueue[here.id] = (c_const_ptrTo(matrix.basis), c_ptrTo(accessor),
                                    c_ptrTo(newRemoteBuffers[0, 0]),
                                    c_ptrTo(newLocalBuffers[0, 0]));
  allLocalesBarrier.barrier();

  const ptrStore : [0 ..# numLocales] globalPtrStoreNoQueue.eltType = globalPtrStoreNoQueue;
  _offDiagInitLocalBuffers(newLocalBuffers, ptrStore);
  _offDiagInitRemoteBuffers(newRemoteBuffers, ptrStore);
  allLocalesBarrier.barrier();
  timers.initialization.stop();

  const numChunks = max(representatives.size * matrix.numberOffDiagTerms() / remoteBufferSize, 1);
  const ranges : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
  const batchedOperatorChunkSize = (representatives.size + numChunks - 1) / numChunks;
  var moreWork : atomic bool = true;
  var curChunkIdx : atomic int = 0;

  logDebug("numChunks = ", numChunks, ", chunkSize = ", batchedOperatorChunkSize);
  coforall taskIdx in 0 ..# numTasks with (ref timers) {
    var batchedOperator = new BatchedOperator(matrix, batchedOperatorChunkSize);

    while moreWork.read() {
      record InTheMeantime {
        inline proc this() {
          return tryProcessLocal(taskIdx, newLocalBuffers, matrix, accessor);
        }
      }
      const inTheMeantime = new InTheMeantime();

      const rangeIdx = curChunkIdx.fetchAdd(1);
      // Multiple threads passed moreWork.read() at once.
      // All whose fetchAdd() was after the one
      // that grabbed the final chunk just break.
      if rangeIdx >= numChunks then
        break;
      // Final rangeIdx -- signal that to everybody
      if rangeIdx == numChunks - 1 then
        moreWork.write(false);
      const r : range(int) = ranges[rangeIdx];

      // Compute a r.size rows of the matrix
      var timer = new Timer();
      timer.start();
      const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
          r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
      var radixOffsets : c_array(int, 257);
      radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
      timer.stop();
      timers.computeOffDiag.add(timer.elapsed(), memoryOrder.relaxed);

      for destLocaleIdx in 0 ..# numLocales {
        const k = radixOffsets[destLocaleIdx];
        const n = radixOffsets[destLocaleIdx + 1] - k;
        ref remoteBuffer = newRemoteBuffers[destLocaleIdx, taskIdx];
        // if n == 0 then continue;
        if destLocaleIdx == here.id {
          assert(!remoteBuffer.isFull.read());
          localProcess(c_const_ptrTo(matrix.basis), c_const_ptrTo(accessor),
                       basisStatesPtr + k,
                       coeffsPtr + k,
                       n);
        }
        else {
          remoteBuffer.submit(basisStatesPtr + k, coeffsPtr + k, n, inTheMeantime);
          inTheMeantime();
        }
      }
    }

    coforall localeIdx in 0 ..# numLocales {
      ref remoteBuffer = newRemoteBuffers[localeIdx, taskIdx];
      while remoteBuffer.isFull.read() {
        if !tryProcessLocal(taskIdx, localeIdx, newLocalBuffers, matrix, accessor) then
          chpl_task_yield();
      }
      remoteBuffer.finish();

      ref localBuffer = newLocalBuffers[localeIdx, taskIdx];
      while !localBuffer.isEOF.read() {
        if !tryProcessLocal(taskIdx, localeIdx, newLocalBuffers, matrix, accessor) then
          chpl_task_yield();
      }
    }
  }

  allLocalesBarrier.barrier();

  for srcLocaleIdx in 0 ..# numLocales {
    for taskIdx in 0 ..# numTasks {
      assert(newLocalBuffers[srcLocaleIdx, taskIdx].isEOF.read());
      assert(newLocalBuffers[srcLocaleIdx, taskIdx].isEmpty.read());
    }
  }
  for srcLocaleIdx in 0 ..# numLocales {
    for taskIdx in 0 ..# numTasks {
      assert(!newRemoteBuffers[srcLocaleIdx, taskIdx].isFull.read());
    }
  }

  timers.put = + reduce newRemoteBuffers._putTime;
  timers.wait = + reduce newRemoteBuffers._waitTime;
  timers.submit = + reduce newRemoteBuffers._submitTime;

  timers.total.stop();
  logDebug("localOffDiagonalNoQueue took ", timers.total.elapsed(), "\n",
           "  ├─ ", timers.initialization.elapsed(), " in initialization\n",
           "  ├─ ", timers.computeOffDiag, " in computeOffDiag\n",
           "  ├─ ", timers.submit, " in remoteBuffer.submit\n",
           "      ├─ ", timers.wait, " in remoteBuffer.isFull.waitFor(false)\n",
           "      └─ ", timers.put, " in remote PUTs\n");
}

private proc localMatrixVector(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                               const ref representatives : [] uint(64)) {
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  localDiagonal(matrix, x, y, representatives);
  if kUseQueue then
    localOffDiagonal(matrix, x, y, representatives);
  else
    localOffDiagonalNoQueue(matrix, x, y, representatives);
}

proc matrixVectorProduct(matrixFilename : string,
                         const ref x,
                         ref y,
                         const ref representatives) {
  coforall loc in Locales with (ref y) do on loc {
    var myMatrix = loadHamiltonianFromYaml(matrixFilename);
    const ref myX = x.getBlock(loc.id)[0, ..];
    ref myY = y.getBlock(loc.id)[0, ..];
    const ref myBasisStates = representatives.getBlock(loc.id);
    myMatrix.basis.uncheckedSetRepresentatives(myBasisStates);
    // myY += 1;
    var timer = new Timer();
    timer.start();
    localMatrixVector(myMatrix, myX, myY, myBasisStates);
    timer.stop();
    logDebug("Spent ", timer.elapsed(), " in localMatrixVector");
  }
}

export proc ls_chpl_matrix_vector_product(matrixPtr : c_ptr(ls_hs_operator), numVectors : c_int,
                                          xPtr : c_ptr(real(64)), yPtr : c_ptr(real(64))) {
  logDebug("Calling ls_chpl_matrix_vector_product ...");
  var matrix = new Operator(matrixPtr, owning=false);
  if matrix.basis.numberWords != 1 then
    halt("bases with more than 64 bits are not yet implemented");
  if matrix.basis.requiresProjection() then
    halt("bases that require projection are not yet supported");
  if numVectors != 1 then
    halt("applying the Operator to more than 1 vector is not yet implemented");

  const ref representatives = matrix.basis.representatives();
  const numStates = representatives.size;
  var x = makeArrayFromPtr(xPtr, (numStates,));
  var y = makeArrayFromPtr(yPtr, (numStates,));
  localMatrixVector(matrix, x, y, representatives);
  logDebug("Done!");
}

/*
proc localMatrixVectorPart(
      H : Operator,
      const ref x : [] ?eltType,
      ref y : [] eltType,
      const ref representatives : [] uint(64)) {
  y = 0;

  // Diagonal contribution
  forall (sigma, i) in zip(representatives, 0..) {

  }


  // Off-diagonal contribution

  coforall loc in Locales do on loc {
    var accessor = new unmanaged ConcurrentAccessor();

    for i,  in zip
    
  }


  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }

  y = 0; // initialize y with zeros
  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }
}
*/

}
