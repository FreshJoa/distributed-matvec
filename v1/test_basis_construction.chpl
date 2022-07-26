use basis;
use Merge;
use MatVec.IO;
use wrapper only initRuntime, deinitRuntime; // , readHDF5Chunk;
use Spawn;

// use BlockDist;
// use BlockCycDist;
// use CyclicDist;
use Distribute;

config const kInputBasisPath = "data/heisenberg_chain_10.yaml";
config const kInputDataPath = "data/heisenberg_chain_10.h5";
config const kOutputDataPath = "output.h5";
config const kChunkSize = 2;
config const kDebugLatticeSymmetries = false;

proc main() {
  initRuntime(kDebugLatticeSymmetries);

  var basis = new DistributedBasis(kInputBasisPath);
  // 1) Construct the representatives in Chapel
  var states = makeStates(basis);
  // 2) Load representatives from HDF5 file and distribute them based on hash
  var expectedStates = 
    distributeStates(loadStates(kInputDataPath, "/basis/representatives"));

  // 1) and 2) should produce the same results
  if !(&& reduce (states.counts == expectedStates.counts)) {
    writeln("[Chapel] ERROR: states.counts and expectedStates.counts do not match!"); 
    return 1;
  }
  for loc in Locales {
    const equal = && reduce (states.representatives[loc.id]
                              == expectedStates.representatives[loc.id]);
    if !equal {
      writeln("[Chapel] ERROR: states.representatives and ",
              "expectedStates.representatives do not match!"); 
      return 1;
    }
  }

  // 3) dump representatives to HDF5 file
  mergeAndWriteStates(states.representatives, states.counts,
                      kOutputDataPath, "/representatives", kChunkSize);
  // the result of 3) should be equal to data in kInputDataPath
  var h5diff = spawn(
    ["h5diff", kOutputDataPath, kInputDataPath, "/representatives", "/basis/representatives"]);
  h5diff.wait();
  if h5diff.exitCode != 0 {
      writeln("[Chapel] mergeAndWriteStates is broken");
      return 1;
  }

  /*
  var N = + reduce states.counts;
  var A = readHDF5Chunk(kInputDataPath, "/basis/representatives", uint(64),
                        (0,), (N,));

  const D : domain(1) dmapped Block(LocaleSpace) = {0 ..# N};
  var A2 : [D] uint(64) = A;

  const x = distributeStates(A2);
  for i in LocaleSpace {
    writeln(x.counts[i], ", ", x.representatives[i]);
  }

  var B = readHDF5Chunk(kInputDataPath, "/hamiltonian/eigenvectors", real(64),
                        (0, 0), (3, N));
  const DB : domain(2) dmapped Block({0 ..# 1, 0 ..# numLocales}) = {0 ..# 3, 0 ..# N};
  var B2 : [DB] real(64) = B;

  writeln(B);
  writeln(B2);

  var C = loadVectors(kInputDataPath, "/hamiltonian/eigenvectors");
  writeln(C[0 ..# 3, ..]);

  const y = distributeVectors(B2, distributionMask(A2));
  for i in LocaleSpace {
    writeln("Chunk #", i);
    writeln(y[i]);
  }
  */

  return 0;
}

proc deinit() {
  deinitRuntime();
}
