use matvec;
use basis;
use states;
use wrapper;
use Distribute;
use Merge;
use MatVec.IO;
use CPtr;
use Time;

config const kInputBasisPath = "data/heisenberg_chain_10.yaml";
config const kInputDataPath = "data/matvec/heisenberg_chain_10.h5";
config const kOutputDataPath = "output.h5";
config const kDebugLatticeSymmetries = false;
config const kPrintOutput = false;
config const kStopEarly = false;
config const kBits = 16;

proc profilingOverview() {
  writeln("╔═════════════════════════╗");
  writeln("║Constructing basis states║");
  writeln("╚═════════════════════════╝");
  _makeStatesTime.print();
  _makeStatesLoopTime.print();
  _makeStatesAllocatingTime.print();
  _makeStatesListToArrayTime.print();
  _makeStatesRemoteCopiesTime.print();
  _splitIntoRangesTime.print();
  writeln("╔═════════════════════╗");
  writeln("║Matrix-vector product║");
  writeln("╚═════════════════════╝");
  _matvecTime.print();
  _batchedApplyTime.print();
  _processBatchTime.print();
  _processBatchInnerTime.print();
  _constructTargetsTime.print();
  // _constructOffsetsTime.print();
  _processLocalTargetsTime.print();
  _processResultsTime.print();
  writeln("╔══════════════╗");
  writeln("║Input / Output║");
  writeln("╚══════════════╝");
  _distributionMaskTime.print();
  _distributionMaskInnerTime.print();
  _kMergeIndicesChunked.print();
  _chunkedIterBodyTime.print();
  _chunkOffsetsTime.print();
  _chunkBoundsTime.print();
  _distributeArrayTime.print();
  _distributeArrayDistributeTime.print();
  _distributeArrayRemoteCopiesTime.print();
  _mergeAndWriteVectorsTime.print();
  _writeHDF5ChunkTime.print();
  _createHDF5DatasetTime.print();
  _loadVectorsTime.print();
  _readHDF5ChunkTime.print();
}

proc main() {
  initRuntime(false);
  var timer = new Timer();
  var basis = new DistributedBasis(kInputBasisPath);
  var states = makeStates(basis, kBits);
  var matrix = new DistributedOperator(kInputBasisPath);
  writeln("[Chapel] Hilbert space dimension is ", states.totalNumberStates());
  writeln("[Chapel] Distribution of states: ", states._counts:real / states.totalNumberStates());

  var mask = distributionMask(states);
  var X = distributeVectors(loadVectors(kInputDataPath, "/x"), mask);

  if (kStopEarly) {
    profilingOverview();
    return 0;
  }

  var Y : [OnePerLocale] [X[0].domain] real(64);
  matvecSimple(matrix, states, X, Y);
  writeln("[Chapel] Done with matvecSimple!");

  var YExpected = distributeVectors(loadVectors(kInputDataPath, "/y"), mask);

  if (kPrintOutput) {
    writeln(Y);
    writeln(YExpected);
  }
  var status : int = 0;
  for loc in Locales {
    const error = max reduce abs(Y[loc.id] - YExpected[loc.id]);
    const norm = max reduce abs(YExpected[loc.id]);
    if (error >= 1e-15 * norm) {
      writeln("[Chapel] ERROR: Y and YExpected do not match!"); 
      status = 1;
    }
    // const equal = && reduce (Y[loc.id] == YExpected[loc.id]);
    // if !equal {
    //   for (y, yExpected) in zip(Y[loc.id], YExpected[loc.id]) {
    //     if (y != yExpected) {
    //       writeln(y, " != ", yExpected, " (Δ = ", yExpected - y, ")");
    //     }
    //   }
    //   writeln("[Chapel] ERROR: Y and YExpected do not match!"); 
    //   status = 1;
    // }
  }

  profilingOverview();
  return status;
}

// proc deinit() {
//   deinitRuntime();
// }
