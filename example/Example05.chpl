use LatticeSymmetries;
use Time;

proc localLoadVectors(filename : string, x : string = "/x", y : string = "/y") {
  var input = readDataset(filename, x, real(64), rank = 2)[0, ..];
  var output = readDataset(filename, y, real(64), rank = 2)[0, ..];
  return (input, output);
}

config const kSystem = "heisenberg_chain";
config const kNumSpins = 10;

proc makeFilename(system, numSpins, extension) {
  return system + "_" + numSpins:string + "." + extension;
}

config const kHamiltonian = "data/" + makeFilename(kSystem, kNumSpins, "yaml");
config const kVectors = if kNumSpins < 24
                          then "data/matvec/" +
                                  makeFilename(kSystem, kNumSpins, "h5")
                          else "data/large-scale/matvec/" +
                                  makeFilename(kSystem, kNumSpins, "h5");
config const kAbsTol = 1e-14;
config const kRelTol = 1e-12;
config const kVerbose = false;

proc approxEqual(a : real, b : real, atol = kAbsTol, rtol = kRelTol) {
  return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)));
}
proc approxEqual(a : [], b : [], atol = kAbsTol, rtol = kRelTol) {
  return [i in a.domain] approxEqual(a[i], b[i], atol, rtol);
}


proc main() {
  initRuntime();
  defer deinitRuntime();

  var timer = new Timer();
  // timer.start();
  var matrix = loadHamiltonianFromYaml(kHamiltonian);
  // timer.stop();
  // logDebug("Reading the Hamiltonian took ", timer.elapsed());

  // timer.clear();
  // timer.start();
  const masks;
  const basisStates = enumerateStates(matrix.basis, masks);
  // timer.stop();
  // logDebug("Enumerating basis states took ", timer.elapsed());

  // for (i, n) in zip(0 .., basisStates._counts) do
  //   writeln("basisStates: ", i, ": ", n, " -> ", n.locale);
  // matrix.basis.build();
  // const (xRaw, yRaw) = localLoadVectors(kVectors);
  // const x = xRaw;
  // const y = xRaw;
  // if kVerbose then
  //   writeln(basisStates);

  // timer.clear();
  // timer.start();
  // const basisStatesRaw = statesFromHashedToBlock(basisStates);
  // timer.stop();
  // logDebug("Merging states took ", timer.elapsed());

  timer.clear();
  timer.start();
  const xBlock = readVectorsAsBlocks(kVectors, "/x");
  timer.stop();
  logDebug("Reading X took ", timer.elapsed());

  // timer.clear();
  // timer.start();
  // const masks = [x in basisStatesRaw] localeIdxOf(x);
  // timer.stop();
  // logDebug("Computing masks took ", timer.elapsed());

  const x = arrFromBlockToHashed(xBlock, masks);
  // timer.clear();
  // timer.start();
  // const x = if kUseNew then arrFromBlockToHashed(masks, xRaw)
  //                      else vectorsFromBlockToHashed(basisStatesRaw, xRaw);
  // timer.stop();
  // logDebug("Distributing X took ", timer.elapsed());

  // for (i, n) in zip(0 .., x._counts) do
  //   writeln("x: ", i, ": ", n, " -> ", n.locale);
  var z = new BlockVector(x.eltType, x.innerDom.dim(0).size, x.counts);
  // for (i, n) in zip(0 .., z._counts) do
  //   writeln("z: ", i, ": ", n, " -> ", n.locale);

  matrixVectorProduct(kHamiltonian, x, z, basisStates);

  const zBlock = arrFromHashedToBlock(z, masks);
  const yBlock = readVectorsAsBlocks(kVectors, "/y");
  // const y = vectorsFromBlockToHashed(basisStatesRaw, yRaw);
  // writeln("y: ", yRaw);
  // writeln("z: ", zRaw);
  const closeEnough = && reduce approxEqual(yBlock, zBlock);
  writeln(closeEnough);
  if (!closeEnough) {
    var maxErrorCount = 10;
    for i in yBlock.domain {
      if !approxEqual(zBlock[i], yBlock[i]) && maxErrorCount > 0 {
        writeln("at ", i, ": ", zBlock[i], " (computed) != ", yBlock[i], " (expected)");
        maxErrorCount -= 1;
      }
    }
  }
  // const x2 = vectorsFromHashedToBlock(basisStates, x);
  // const y = readVectorsAsBlocks(kVectors, "/y");
  // writeln(x);

// proc matrixVectorProduct(matrixFilename : string,
//                          const ref x,
//                          ref y,
//                          const ref representatives) {

  // var timer = new Timer();
  // timer.start();
  // localMatVec(matrix, x, z);
  // timer.stop();

  // const closeEnough = && reduce approxEqual(y, z);
  // if (!closeEnough) {
  //   var maxErrorCount = 10;
  //   for i in x.domain {
  //     if !approxEqual(z[i], y[i]) && maxErrorCount > 0 {
  //       writeln("at ", i, ": ", z[i], " != ", y[i]);
  //       maxErrorCount -= 1;
  //     }
  //   }
  // }
  // assert(closeEnough);
  // writeln(timer.elapsed());

  // const reference = readBasisStatesAsBlocks(kRepresentatives, "/basis/representatives");
  // for loc in Locales {
  //   const i = reference.domain.localSubdomain(loc);
  //   writeln(loc, ": ", reference[i]);
  // }

  // const basisStates = statesFromBlockToHashed(reference);
  // assert(reference.size == + reduce basisStates._counts);
  // writeln(basisStates);

  // const obtained = statesFromHashedToBlock(basisStates);
  // assert(obtained.size == reference.size);
  // writeln(obtained == reference);
  // writeln(&& reduce (obtained == reference));

  // for loc in Locales {
  //   const i = obtained.domain.localSubdomain(loc);
  //   writeln(loc, ": ", obtained[i]);
  // }
}
