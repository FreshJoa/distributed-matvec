use LatticeSymmetries;
use Time;

proc localLoadVectors(filename : string, x : string = "/x", y : string = "/y") {
  var input = readDataset(filename, x, real(64), rank = 2)[0, ..];
  var output = readDataset(filename, y, real(64), rank = 2)[0, ..];
  return (input, output);
}

config const kHamiltonian = "data/heisenberg_chain_10.yaml";
config const kVectors = "data/matvec/heisenberg_chain_10.h5";
config const kAbsTol = 1e-14;
config const kRelTol = 1e-12;

proc approxEqual(a : real, b : real, atol = kAbsTol, rtol = kRelTol) {
  return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)));
}
proc approxEqual(a : [], b : [], atol = kAbsTol, rtol = kRelTol) {
  return [i in a.domain] approxEqual(a[i], b[i], atol, rtol);
}

proc main() {
  ls_hs_init();

  var matrix = loadHamiltonianFromYaml(kHamiltonian);
  matrix.basis.build();
  const (x, y) = localLoadVectors(kVectors);
  var z : x.type;

  var timer = new Timer();
  timer.start();
  localMatVec(matrix, x, z);
  timer.stop();

  const closeEnough = && reduce approxEqual(y, z);
  if (!closeEnough) {
    var maxErrorCount = 10;
    for i in x.domain {
      if !approxEqual(z[i], y[i]) && maxErrorCount > 0 {
        writeln("at ", i, ": ", z[i], " != ", y[i]);
        maxErrorCount -= 1;
      }
    }
  }
  assert(closeEnough);
  writeln(timer.elapsed());
}

proc deinit() {
  ls_hs_exit();
}
