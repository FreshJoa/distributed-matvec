use basis;
use Distribute;
use MatVec.IO;
use wrapper only initRuntime, deinitRuntime;
use Spawn;

config const kInputDataPath = "/home/tom/src/spin-ed/example/data/heisenberg_kagome_12.h5";
config const kOutputDataPath = "output.h5";
config const kDebugLatticeSymmetries = false;

proc main() {
  initRuntime(kDebugLatticeSymmetries);

  const blockDistRepresentatives = loadStates(kInputDataPath, "/basis/representatives");
  const blockDistVectors = loadVectors(kInputDataPath, "/hamiltonian/eigenvectors", real(64));
  const mask = distributionMask(blockDistRepresentatives);

  var states = distributeStates(blockDistRepresentatives, mask);
  var vectors = distributeVectors(blockDistVectors, mask);

  mergeAndWriteVectors(states.representatives, vectors, states.counts, kOutputDataPath, "/y");
  var h5diff = spawn(
    ["h5diff", kOutputDataPath, kInputDataPath, "/y", "/hamiltonian/eigenvectors"]);
  h5diff.wait();
  if h5diff.exitCode != 0 {
      writeln("[Chapel] mergeAndWriteVectors is broken");
      return 1;
  }
  return 0;
}

proc deinit() {
  deinitRuntime();
}

