use LatticeSymmetries;
use Time;

proc localLoadRepresentatives(filename : string, dataset : string = "/basis/representatives") {
  return readDataset(filename, dataset, uint(64), rank = 1);
}

config const kBasis = "data/heisenberg_kagome_12.yaml";
config const kRepresentatives = "data/heisenberg_kagome_12.h5";

proc main() {
  ls_hs_init();

  const reference = localLoadRepresentatives(kRepresentatives);
  var basis = loadBasisFromYaml(kBasis);
  var timer = new Timer();
  timer.start();
  basis.build();
  timer.stop();
  const ref predicted = basis.representatives();
  const theSame = && reduce [i in reference.domain] reference[i] == predicted[i];
  if !theSame {
    ref predicted = basis.representatives();
    var maxErrorCount = 10;
    for i in reference.domain {
      if reference[i] != predicted[i] && maxErrorCount > 0 {
        writeln("at index ", i, ": ", reference[i], " != ", predicted[i]);
        maxErrorCount -= 1;
      }
    }
  }
  writeln(timer.elapsed());
}

proc deinit() {
  ls_hs_exit();
}
