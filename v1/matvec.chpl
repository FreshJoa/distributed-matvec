
use basis;
use wrapper;
use Search;
use CPtr;
use Time;
use states only hash64_01;
use profiling;
use DynamicIters;
use CommDiagnostics;
use Memory.Diagnostics;

config const kBatchSize = 100;

var _batchedApplyTime = new MeasurementTable("DistributedOperator::batchedApply");

class DistributedOperator {
  var _localBases : [OnePerLocale] ls_hs_spin_basis_v1;
  var _localOperators : [OnePerLocale] ls_hs_operator_v1;

  proc init(path: string) {
    this._localBases = [loc in Locales] new ls_hs_spin_basis_v1(nil, nil);
    this._localOperators = [loc in Locales] new ls_hs_operator_v1(nil, nil);
    complete();

    // for o in this._localOperators {
    //   writeln(o);
    // }
    forall loc in Locales do on loc {
      ls_hs_basis_and_hamiltonian_from_yaml(
        path.localize().c_str(),
        c_ptrTo(this._localBases[loc.id]),
        c_ptrTo(this._localOperators[loc.id]));
    }
    // for o in this._localOperators {
    //   writeln(o);
    // }
  }

  proc deinit() {
    forall loc in Locales do on loc {
      // writeln("ls_hs_destroy_operator on ", this._localOperators[loc.id]);
      ls_hs_destroy_operator(c_ptrTo(this._localOperators[loc.id]));
      // writeln("ls_hs_destroy_basis on ", this._localBases[loc.id]);
      ls_hs_destroy_spin_basis(c_ptrTo(this._localBases[loc.id]));
      // writeln("done: operator=", this._localOperators[loc.id],
      //         " basis=", this._localBases[loc.id]);
    }
    // writeln("[Chapel] DistributedOperator.deinit is done");
  }

  proc rawPtr() { return this._localOperators[here.id].payload; }
  proc maxEntriesPerRow() { return ls_operator_max_buffer_size(rawPtr()):int; }
  proc requiredBufferSize(count : int = 1) {
    return count * maxEntriesPerRow();
  }

  proc batchedApply(
    spins : [?D] uint(64),
    outOffsets : [] uint(64),
    outSpins : [] uint(64),
    outCoeffs : [] complex(128),
    offset : int = 0,
    count : int = spins.size)
     where spins.domain.rank == 1 &&
           outOffsets.domain.rank == 1 &&
           outSpins.domain.rank == 1 &&
           outCoeffs.domain.rank == 1 {
    var __timer = getTimerFor(_batchedApplyTime);
    const status = ls_hs_operator_apply(
      c_ptrTo(this._localOperators[here.id]),
      count:uint(64),
      c_ptrTo(spins[spins.domain.low + offset]),
      c_ptrTo(outOffsets),
      c_ptrTo(outSpins),
      c_ptrTo(outCoeffs));
    assert(status == 0);
  }
}

record TaskState {
  const batchSize : int;
  const bufferSize : int;
  var offsetsBatch : [0 ..# batchSize + 1] uint(64);
  var spinsBatch : [0 ..# bufferSize] uint(64);
  var coeffsBatch : [0 ..# bufferSize] complex(128);
  var ysBatch : [0 ..# batchSize] complex(128);
}

// var processBatchTimer = new MeasurementTable("processBatch");
// var processBatchTimer_1 = new MeasurementTable("processBatch::1");
// var processBatchTimer_2 = new MeasurementTable("processBatch::2");
// var processBatchTimer_3 = new MeasurementTable("processBatch::3");
// var processBatchTimer_4 = new MeasurementTable("processBatch::4");
// var copyToHere = new MeasurementTable("copyToHere");
// var getIndexTimer = new MeasurementTable("getIndex");

var _processBatchTime = new MeasurementTable("BatchProcessor::process");
var _processBatchInnerTime = new MeasurementTable("BatchProcessor::process::inner loop");
var _constructTargetsTime = new MeasurementTable("BatchProcessor::_constructTargetsTimer");
var _constructOffsetsTime = new MeasurementTable("BatchProcessor::_constructOffsetsTimer");
var _processLocalTargetsTime = new MeasurementTable("BatchProcessor::_processLocalTargetsTimer");
var _processResultsTime = new MeasurementTable("BatchProcessor::_processResultsTimer");

record BatchProcessor {
  var bufferSize : int;
  var targets : [LocaleSpace] [0 ..# bufferSize] (int, uint(64), complex(128));
  var targetSizes : [LocaleSpace] int;
  var offsets : [0 ..# numLocales + 1] int;
  var results : [0 ..# bufferSize] (int, complex(128));

  proc init(bufferSize : int) {
    this.bufferSize = bufferSize;
  }

  proc _constructTargets(count : int, const ref info : TaskState) {
    var __timer = getTimerFor(_constructTargetsTime);
    targetSizes = 0;
    for k in 0 ..# count {
      for j in info.offsetsBatch[k]:int .. info.offsetsBatch[k + 1]:int - 1 {
        const sj = info.spinsBatch[j];
        const cj = info.coeffsBatch[j];
        const localeId = (hash64_01(sj) % numLocales:uint):int;
        targets[localeId][targetSizes[localeId]] = (k, sj, cj);
        targetSizes[localeId] += 1;
      }
    }
  }

  proc _constructOffsets() {
    assert(offsets.locale == here);
    assert(targetSizes.locale == here);
    var n : int = 0;
    offsets[0] = n;
    for locId in LocaleSpace {
      n += targetSizes[locId];
      offsets[locId + 1] = n;
    }
  }

  proc _processLocalTargets(const ref localTargets : [] (int, uint(64), complex(128)),
                            ref localResults : [] (int, complex(128)),
                            const ref basisStates,
                            const ref X) {
    var __timer = getTimerFor(_processLocalTargetsTime);
    const ref localX = X[here.id];
    const getterFn = basisStates.getIndexFn();
    // startVerboseComm();
    foreach ((k, sj, cj), r) in zip(localTargets, localResults) {
      const j = getterFn(sj); // basisStates.getIndex(sj);
      const xj = localX[0, j];
      r = (k, cj * xj);
    }
    // stopVerboseComm();
  }

  proc _processResults(ref info : TaskState) {
    var __timer = getTimerFor(_processResultsTime);
    info.ysBatch = 0;
    foreach (k, yk) in results[0 ..# offsets[numLocales]] {
      info.ysBatch[k] += yk;
    }
  }

  proc process(count : int,
               ref info : TaskState,
               const ref basisStates : BasisStates,
               const ref X) {
    var __timer = getTimerFor(_processBatchTime);
    _constructTargets(count, info);
    _constructOffsets();
    var __timer1 = getTimerFor(_processBatchInnerTime);
    coforall loc in Locales with (ref _processLocalTargetsTime) do on loc {
      const localTargets = targets[loc.id][0 ..# targetSizes[loc.id]];
      const localOffsets = offsets;
      // startVerboseComm();
      var localResults : [0 ..# localTargets.size] (int, complex(128));
      _processLocalTargets(localTargets, localResults, basisStates, X);
      // stopVerboseComm();
      // assert(offsets[loc.id + 1] - offsets[loc.id] == localResults.size);
      // assert(offsets[loc.id + 1] <= results.size);
      results[localOffsets[loc.id] .. localOffsets[loc.id + 1] - 1] = localResults;
    }
    stopVerboseComm();
    __timer1.stop();
    _processResults(info);
  }
};

// proc processBatch(count : int,
//                   ref info : TaskState,
//                   const ref basisStates : BasisStates,
//                   const ref X) {
//   var _timer = new ScopedTimer(processBatchTimer);
//   type eltType = X[0].eltType;
// 
//   var _timer_1 = new ScopedTimer(processBatchTimer_1);
//   var targets : [0 ..# numLocales] [0 ..# info.bufferSize] (int, uint(64));
//   var targetsSizes : [LocaleSpace] int;
//   var results : [0 ..# info.bufferSize] real;
//   _timer_1.stop();
//   var _timer_2 = new ScopedTimer(processBatchTimer_2);
//   for j in 0 ..# info.offsetsBatch[count]:int {
//     const sj = info.spinsBatch[j];
//     const xj = info.coeffsBatch[j];
//     const localeId = (hash64_01(sj) % numLocales:uint):int;
//     targets[localeId][targetsSizes[localeId]] = (j, sj);
//     targetsSizes[localeId] += 1;
//   }
//   _timer_2.stop();
//   var _timer_3 = new ScopedTimer(processBatchTimer_3);
//   // startVerboseComm();
//   coforall loc in Locales with (ref copyToHere, ref getIndexTimer) do on loc {
//     var _t_1 = new ScopedTimer(copyToHere);
//     const localTargets = targets[loc.id][0 ..# targetsSizes[loc.id]];
//     _t_1.stop();
//     var _t_2 = new ScopedTimer(getIndexTimer);
//     for (i, sj) in localTargets {
//       results[i] = X[loc.id][0, basisStates.getIndex(sj)];
//     }
//   }
//   // stopVerboseComm();
//   _timer_3.stop();
// 
//   var _timer_4 = new ScopedTimer(processBatchTimer_4);
//   for k in 0 ..# count {
//     var yk : complex(128) = 0;
//     for _j in info.offsetsBatch[k] .. info.offsetsBatch[k + 1] - 1 {
//       const xj = results[_j:int];
//       const cj = info.coeffsBatch[_j:int];
//       yk += cj * xj;
//     }
//     info.ysBatch[k] = yk;
//   }
//   _timer_4.stop();
// }
// 
// inline proc batchedApply(matrix : DistributedOperator, ref spins, ref s : TaskState) {
//   matrix.batchedApply(
//     spins,
//     s.offsetsBatch,
//     s.spinsBatch,
//     s.coeffsBatch);
// }
var _matvecTime = new MeasurementTable("matvecSimple");

proc matvecSimple(matrix : DistributedOperator, 
                  const ref basisStates : BasisStates,
                  const ref X, // [OnePerLocale] [?D] ?eltType
                  ref Y) // [OnePerLocale] [?D] ?eltType
{
  var __timer = getTimerFor(_matvecTime);

  const batchSize = kBatchSize;
  const bufferSize = matrix.requiredBufferSize(batchSize);
  type eltType = X[0].eltType;
  writeln("[Chapel] batchSize = ", batchSize, ", bufferSize = ", bufferSize);
  // startVerboseMem();
  // startVerboseComm();
  coforall loc in Locales do on loc {
    ref representatives = basisStates.getRepresentatives()[0 ..# basisStates.getCounts()];
    const numberChunks = (representatives.size + batchSize - 1) / batchSize;
    writeln("[Chapel] locale ", loc.id, " has ", numberChunks, " chunks to process");
    var info = new TaskState(batchSize, bufferSize);
    var batchProcessor = new BatchProcessor(info.bufferSize);

    forall chunkId in dynamic(0 ..# numberChunks) with (in info, in batchProcessor) {
      const i = chunkId * batchSize;
      const n = min(batchSize, representatives.size - i);
      matrix.batchedApply(
        representatives,
        info.offsetsBatch,
        info.spinsBatch,
        info.coeffsBatch,
        /*offset=*/i,
        /*count=*/n);
      batchProcessor.process(n, info, basisStates, X);

      if isSubtype(Y[loc.id].eltType, complex) {
        Y[loc.id][0, i ..# n] = info.ysBatch[0 ..# n];
      }
      else {
        assert(&& reduce (info.ysBatch[0 ..# n].im == 0));
        Y[loc.id][0, i ..# n] = info.ysBatch[0 ..# n].re;
        // for k in 0 ..# n {
        //   assert(info.ysBatch[k].im == 0);
        //   Y[loc.id][0, i + k] = info.ysBatch[k].re;
        // }
      }
    }
  }
  // stopVerboseComm();
  // stopVerboseMem();

}













