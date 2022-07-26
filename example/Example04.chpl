use LatticeSymmetries;
use Time;

config const kBasis = "data/heisenberg_kagome_12.yaml";
config const kRepresentatives = "data/construction/heisenberg_kagome_12.h5";

proc main() {
  initRuntime();
  defer deinitRuntime();

  const reference = readBasisStatesAsBlocks(kRepresentatives, "/basis/representatives");
  // for loc in Locales {
  //   const i = reference.domain.localSubdomain(loc);
  //   writeln(loc, ": ", reference[i]);
  // }

  const basisStates = statesFromBlockToHashed(reference);
  assert(reference.size == + reduce basisStates._counts);
  // writeln(basisStates);

  const obtained = statesFromHashedToBlock(basisStates);
  assert(obtained.size == reference.size);
  // writeln(obtained == reference);
  writeln(&& reduce (obtained == reference));

  // for loc in Locales {
  //   const i = obtained.domain.localSubdomain(loc);
  //   writeln(loc, ": ", obtained[i]);
  // }
}
