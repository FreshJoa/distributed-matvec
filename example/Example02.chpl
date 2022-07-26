use LatticeSymmetries;
use Time;
use RangeChunk;

config const kNumSpins = 10;
config const kVerbose = false;

proc main() {
  initRuntime();
  defer deinitRuntime();

  var basis0 = new Basis("{\"number_spins\": " + kNumSpins:string
                        + ", \"hamming_weight\": " + (kNumSpins / 2):string
                        + "}");
  const masks0;
  const buckets0 = enumerateStates(basis0, masks0);

  const batchSize = 3;
  var fakeVectors = new BlockVector(real(64), batchSize, buckets0.counts);
  coforall loc in Locales with (ref fakeVectors) do on loc {
    const ref src = buckets0[loc];
    ref dest = fakeVectors[loc];
    for batchIdx in 0 ..# batchSize {
      dest[batchIdx, ..] = src;
    }
  }

  if kVerbose then
    for loc in Locales do
      on loc do
        logDebug(buckets0[loc]);
  logDebug("total number of states: ", + reduce buckets0.counts);
  if kVerbose then
    logDebug("masks: ", masks0);

  var block = arrFromHashedToBlock(fakeVectors, masks0);
  if kVerbose then
    for loc in Locales do
      on loc do
        logDebug(block[block.localSubdomain()]);

  var fakeVectors2 = arrFromBlockToHashed(block, masks0);
  if kVerbose then
    for loc in Locales do
      on loc do
        logDebug(fakeVectors2[loc]);

  return 0;
}
