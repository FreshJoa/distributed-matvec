module StatesEnumeration {

use FFI;
use ForeignTypes;
use Vector;

use BitOps;
use BlockDist;
use ChplConfig;
use CTypes;
use CyclicDist;
use DynamicIters;
use IO;
use RangeChunk;
use Search;
use AllLocalesBarriers;
use CommDiagnostics;
use Time;
use Memory.Initialization;

config const isRepresentativeBatchSize : int = 65536;
config const enumerateStatesNumChunks : int = 2 * numLocales * here.maxTaskPar;
config const statesFromHashedToBlockNumChunks = 2 * here.maxTaskPar;
config const statesFromBlockToHashedNumChunks = 2 * numLocales * here.maxTaskPar;
config const kUseLowLevelComm : bool = true;
// config const numChunksPerLocale = 3;

inline proc GET(addr, node, rAddr, size) {
  __primitive("chpl_comm_get", addr, node, rAddr, size);
}

inline proc PUT(addr, node, rAddr, size) {
  __primitive("chpl_comm_put", addr, node, rAddr, size);
}

/* Get the next integer with the same Hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
private inline proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

private inline proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

private inline proc nextState(v: uint(64), param isHammingWeightFixed : bool): uint(64) {
  return if isHammingWeightFixed then nextStateFixedHamming(v) else nextStateGeneral(v);
}

private inline iter _manyNextStateIter(in v: uint(64), bound: uint(64),
                                       param isHammingWeightFixed : bool) {
  while true {
    yield v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if v == bound then break;
    v = nextState(v, isHammingWeightFixed);
  }
}
private inline iter manyNextState(v: uint(64), bound: uint(64),
                                  isHammingWeightFixed : bool) {
  assert(v <= bound, "v is greater than bound");
  if isHammingWeightFixed then
    for x in _manyNextStateIter(v, bound, true) do
      yield x;
  else
    for x in _manyNextStateIter(v, bound, false) do
      yield x;
}
private proc manyNextState(v: uint(64), bound: uint(64), ref buffer: [] uint(64),
                           isHammingWeightFixed : bool) : int {
  if buffer.size == 0 then return 0;

  var offset = 0;
  for (x, i) in zip(manyNextState(v, bound, isHammingWeightFixed), 0 ..) {
    buffer[offset] = x;
    offset += 1;
    if offset == buffer.size then break;
  }

  return offset;
}


/*
private inline iter findStatesInRange(in lower: uint(64), upper: uint(64),
                                      isHammingWeightFixed : bool,
                                      const ref basis : Basis) {
  if lower > upper then
    return;
  var buffer: [0 ..# isRepresentativeBatchSize] uint(64) = noinit;
  var flags: [0 ..# isRepresentativeBatchSize] uint(8) = noinit;
  var norms: [0 ..# isRepresentativeBatchSize] real(64) = noinit;
  while (true) {
    const written = manyNextState(lower, upper, buffer, isHammingWeightFixed);
    isRepresentative(basis, buffer[0 ..# written], flags[0 ..# written], norms[0 ..# written]);
    for i in 0 ..# written {
      if (flags[i]:bool && norms[i] > 0) { yield buffer[i]; }
    }
    if (buffer[written - 1] == upper) { break; }

    if isHammingWeightFixed { lower = nextState(buffer[written - 1], true); }
    else { lower = nextState(buffer[written - 1], false); }
  }
}
*/

private inline proc unprojectedStateToIndex(basisState : uint(64), isHammingWeightFixed : bool) : int {
  return if isHammingWeightFixed then ls_hs_fixed_hamming_state_to_index(basisState)
                                 else basisState:int;
}
private inline proc unprojectedIndexToState(stateIndex : int, hammingWeight : int) : uint(64) {
  return if hammingWeight != -1 then ls_hs_fixed_hamming_index_to_state(stateIndex, hammingWeight:c_int)
                                else stateIndex:uint(64);
}

/* We want to split `r` into `numChunks` non-overlapping ranges but such that for each range
   `chunk` we have that `popcount(chunk.low) == popcount(chunk.high) == popcount(r.low)` if
   `isHammingWeightFixed` was set to `true`.
 */
proc determineEnumerationRanges(r : range(uint(64)), in numChunks : int,
                                isHammingWeightFixed : bool) {
  var timer = new Timer();
  timer.start();
  const hammingWeight = if isHammingWeightFixed then popcount(r.low):int else -1;
  const lowIdx = unprojectedStateToIndex(r.low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(r.high, isHammingWeightFixed);
  const totalSize = highIdx - lowIdx + 1;
  if numChunks > totalSize then
    numChunks = totalSize;
  var ranges : [0 ..# numChunks] range(uint(64));
  // var ranges = newCyclicArr({0 ..# numChunks}, range(uint(64)));
  for (r, i) in zip(chunks(lowIdx .. highIdx, numChunks), 0..) {
    ranges[i] = unprojectedIndexToState(r.low, hammingWeight)
                  .. unprojectedIndexToState(r.high, hammingWeight);
  }
  timer.stop();
  // logDebug("determineEnumerationRanges(", r, ", ", numChunks, ", ", isHammingWeightFixed,
  //          ") took ", timer.elapsed());
  return ranges;
}

/* Hash function which we use to map spin configurations to locale indices.

   Typical usage:
   ```chapel
   var localeIndex = (hash64_01(x) % numLocales:uint):int;
   ```
*/
private inline proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

inline proc localeIdxOf(basisState : uint(64)) : int 
    where CHPL_COMM == "" {
  return 0;
}
inline proc localeIdxOf(basisState : uint(64)) : int 
    where CHPL_COMM != "" {
  return (hash64_01(basisState) % numLocales:uint):int;
}

private proc _enumStatesComputeMasksAndCounts(const ref states, ref outMasks) {
  const totalCount = states.size;
  var counts : [0 ..# numLocales] int;

  outMasks.resize(totalCount);
  foreach i in 0 ..# totalCount {
    const key = localeIdxOf(states[i]);
    outMasks[i] = key:uint(8);
    counts[key] += 1;
  }
  return counts;
}

private proc _enumerateStatesProjected(r : range(uint(64)), const ref basis : Basis,
                                       ref outStates : Vector(uint(64))) {
  if r.size == 0 then return;
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  var buffer: [0 ..# isRepresentativeBatchSize] uint(64) = noinit;
  var flags: [0 ..# isRepresentativeBatchSize] uint(8) = noinit;
  var norms: [0 ..# isRepresentativeBatchSize] real(64) = noinit;
  var lower = r.low;
  const upper = r.high;
  while true {
    const written = manyNextState(lower, upper, buffer, isHammingWeightFixed);
    isRepresentative(basis, buffer[0 ..# written],
                            flags[0 ..# written],
                            norms[0 ..# written]);
    for i in 0 ..# written do
      if flags[i]:bool && norms[i] > 0 then
        outStates.pushBack(buffer[i]);

    const last = buffer[written - 1];
    if last == upper then break;

    if isHammingWeightFixed { lower = nextState(last, true); }
    else { lower = nextState(last, false); }
  }
}
private proc _enumerateStatesUnprojected(r : range(uint(64)), const ref basis : Basis,
                                         ref outStates : Vector(uint(64))) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  if isHammingWeightFixed && popcount(r.low) != popcount(r.high) then
    halt("r.low=" + r.low:string + " and r.high=" + r.high:string
        + " have different Hamming weight: " + popcount(r.low):string
        + " vs. " + popcount(r.high):string);
  const lowIdx = unprojectedStateToIndex(r.low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(r.high, isHammingWeightFixed);
  const totalCount = highIdx - lowIdx + 1;
  outStates.resize(totalCount);
  manyNextState(r.low, r.high, outStates._arr, isHammingWeightFixed);
  // return _enumStatesComputeMasksAndCounts(outStates, outMasks);
}
private proc _enumerateStatesUnprojectedSpinfulFermion(r : range(uint(64)),
                                                       const ref basis : Basis,
                                                       ref outStates : Vector(uint(64))) {
  assert(basis.isSpinfulFermionicBasis());
  const numberSites = basis.numberSites();
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basis.numberUp();
  const numberDown = basis.numberParticles() - numberUp;

  const minA = r.low & mask;
  const maxA = r.high & mask;
  const countA =
    unprojectedStateToIndex(maxA, isHammingWeightFixed=true) -
      unprojectedStateToIndex(maxA, isHammingWeightFixed=true) + 1;

  const minB = (r.low >> numberSites) & mask;
  const maxB = (r.high >> numberSites) & mask;
  const countB =
    unprojectedStateToIndex(maxB, isHammingWeightFixed=true) -
      unprojectedStateToIndex(minB, isHammingWeightFixed=true) + 1;

  outStates.resize(countA * countB);
  var offset = 0;
  for xB in manyNextState(minB, maxB, isHammingWeightFixed=true) {
    for xA in manyNextState(minA, maxA, isHammingWeightFixed=true) {
      const x = (xB << numberSites) | xA;
      outStates[offset] = x;
      offset += 1;
    }
  }
}
private proc _enumerateStates(r : range(uint(64)), const ref basis : Basis,
                              ref outStates : Vector(uint(64))) {
  if basis.requiresProjection() then
    _enumerateStatesProjected(r, basis, outStates);
  else if basis.isStateIndexIdentity() || basis.isHammingWeightFixed() then
    _enumerateStatesUnprojected(r, basis, outStates);
  else
    _enumerateStatesUnprojectedSpinfulFermion(r, basis, outStates);
}

/*
private proc _enumerateStatesUnprojected(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  if isHammingWeightFixed && popcount(r.low) != popcount(r.high) then
    halt("r.low=" + r.low:string + " and r.high=" + r.high:string
        + " have different Hamming weight: " + popcount(r.low):string
        + " vs. " + popcount(r.high):string);
  const lowIdx = unprojectedStateToIndex(r.low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(r.high, isHammingWeightFixed);
  const totalCount = highIdx - lowIdx + 1;

  const estimatedCountPerLocale = totalCount / numLocales;
  foreach i in LocaleSpace do
    outVectors[i].reserve(estimatedCountPerLocale);

  var v = r.low;
  for i in 0 ..# totalCount {
    const localeIdx = localeIdxOf(v);
    outVectors[localeIdx].pushBack(v);
    v = if isHammingWeightFixed
          then nextStateFixedHamming(v)
          else nextStateGeneral(v);
  }
}
*/
/*
private proc _enumerateStatesSpinfulFermion(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  assert(basis.isSpinfulFermionicBasis());
  const numberSites = basis.numberSites();
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basis.numberUp();
  const numberDown = basis.numberParticles() - numberUp;
  const countA =
    unprojectedStateToIndex(r.high & mask, isHammingWeightFixed=true) -
      unprojectedStateToIndex(r.low & mask, isHammingWeightFixed=true) + 1;
  const countB =
    unprojectedStateToIndex((r.high >> numberSites) & mask, isHammingWeightFixed=true) -
      unprojectedStateToIndex((r.low >> numberSites) & mask, isHammingWeightFixed=true) + 1;
  const estimatedCountPerLocale = countA * countB / numLocales;
  foreach i in LocaleSpace do
    outVectors[i].reserve(estimatedCountPerLocale);

  var vB = (r.low >> numberSites) & mask;
  for _iB in 0 ..# countB {
    var vA = r.low & mask;
    for _iA in 0 ..# countA {
      const v = (vB << numberSites) | vA;
      const localeIdx = localeIdxOf(v);
      outVectors[localeIdx].pushBack(v);

      vA = nextStateFixedHamming(vA);
    }
    vB = nextStateFixedHamming(vB);
  }
}
*/
/*
private proc _enumerateStatesProjected(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  for v in findStatesInRange(r.low, r.high, isHammingWeightFixed, basis) {
    const localeIdx = localeIdxOf(v);
    outVectors[localeIdx].pushBack(v);
  }
}
*/
/*
private proc _enumerateStates(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  // logDebug("_enumerateStates");
  if basis.requiresProjection() {
    _enumerateStatesProjected(r, basis, outVectors);
  }
  else if basis.isStateIndexIdentity() || basis.isHammingWeightFixed() {
    _enumerateStatesUnprojected(r, basis, outVectors);
  }
  else {
    _enumerateStatesSpinfulFermion(r, basis, outVectors);
  }
}
*/

/*
proc _emptyVectors(type t)
  // : [LocaleSpace] Vector(t)
  { return [loc in LocaleSpace] new Vector(t); }


proc enumerateStates(const ref basis : Basis, in numChunks : int, globalRange : range(uint(64))) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  const ranges = determineEnumerationRanges(globalRange, numChunks, isHammingWeightFixed);
  assert(ranges.size <= numChunks);
  numChunks = ranges.size; // We could have fewer elements in the range than numChunks.
                           // In such case, we limit numChunks

  const D = {0 ..# numChunks} dmapped Cyclic(startIdx=0);
  var buckets: [D] [LocaleSpace] Vector(uint(64)) = [d in D] _emptyVectors(uint(64));
  // noinit;
  // new Vector(uint(64));

  startVerboseComm();
  forall i in D {
    // (r, i) in zip(ranges, 0..) {
    const r = ranges[i];
    assert(here == buckets[i].locale);
    const localBasis = basis;
    _enumerateStates(r, localBasis, buckets[i]);
  }
  stopVerboseComm();

  const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
  var representatives : [OnePerLocale] Vector(uint(64)) =
    [loc in OnePerLocale] new Vector(uint(64));
  forall (vector, localeIdx) in zip(representatives, LocaleSpace) {
    vector.reserve(+ reduce [bucket in buckets] bucket[localeIdx].size);
    for bucket in buckets {
      const ref v = bucket[localeIdx];
      vector.append(v);
    }
  }
  return new BlockVector(representatives);
}
proc enumerateStates(const ref basis : Basis, numChunks : int = enumerateStatesNumChunks) {
  const lower = basis.minStateEstimate();
  const upper = basis.maxStateEstimate();
  return enumerateStates(basis, numChunks, lower .. upper);
}
*/

proc prefixSum(arr : [] ?eltType) {
  var sums : [0 ..# arr.size] eltType;
  var total : eltType = 0;
  for i in sums.domain {
    sums[i] = total;
    total += arr[i];
  }
  return sums;
}

proc permuteBasedOnMasks(arrSize : int, masks : c_ptr(?maskType), arr : c_ptr(?eltType),
                         counts : [] int, destOffsets : [] int,
                         destPtrs : [] c_ptr(eltType)) {
  var offsets : [0 ..# numLocales] int = prefixSum(counts);
  var src : [0 ..# arrSize] eltType = noinit;
  for i in 0 ..# arrSize {
    const key = masks[i]:int;
    src[offsets[key]] = arr[i];
    offsets[key] += 1;
  }

  var copyTimer = new Timer();
  copyTimer.start();
  var i = 0;
  for localeIdx in 0 ..# numLocales {
    if counts[localeIdx] > 0 {
      const srcPtr = c_ptrTo(src[i]);
      const destPtr = destPtrs[localeIdx] + destOffsets[localeIdx];
      const destSize = counts[localeIdx]:c_size_t * c_sizeof(eltType);
      PUT(srcPtr, localeIdx, destPtr, destSize);
      i += counts[localeIdx];
    }
  }
  copyTimer.stop();

  return copyTimer.elapsed();
}

inline proc _enumStatesMakeBuckets(numChunks : int) {
  const dom = {0 ..# numChunks} dmapped Cyclic(startIdx=0);
  var buckets : [dom] (Vector(uint(64)), Vector(uint(8)));
  return buckets;
}

proc _enumStatesComputeCounts(ref buckets,
                              const ref ranges : [] range(uint(64)),
                              const ref basis : Basis) {
  assert(!basis.requiresProjection());
  assert(basis.isStateIndexIdentity() || basis.isHammingWeightFixed());
  assert(here.id == 0);

  const numChunks = ranges.size;
  var counts : [0 ..# numChunks, 0 ..# numLocales] int;
  const countsPtr = c_ptrTo(counts[0, 0]);

  // if kVerboseComm then startVerboseComm();
  // var countsTimer = new Timer();
  // countsTimer.start();
  const serializedBasis : string = basis.toJSON();
  coforall loc in Locales do on loc {
    const myBasis = new Basis(serializedBasis.localize());
    const myRanges : [0 ..# numChunks] range(uint(64)) = ranges;
    const mySubdomain = buckets.localSubdomain();

    forall chunkIdx in mySubdomain {
      ref (outStates, outMasks) = buckets.localAccess(chunkIdx);
      // This is the actual computation!
      _enumerateStates(myRanges[chunkIdx], myBasis, outStates);
      const myCounts = _enumStatesComputeMasksAndCounts(outStates, outMasks);
      // const myCounts = _enumerateStatesUnprojectedNew(myRanges[chunkIdx], myBasis,
      //                                                 outStates, outMasks);
      // Copy computed chunks back to locale 0. We use low-level PUT, because
      // Chapel generates bad code otherwise
      const myCountsPtr = c_const_ptrTo(myCounts[0]);
      const destPtr = countsPtr + chunkIdx * numLocales;
      const destSize = numLocales:c_size_t * c_sizeof(int);
      PUT(myCountsPtr, 0, destPtr, destSize);
    }
  }
  // countsTimer.stop();
  // if kVerboseComm then stopVerboseComm();

  return counts;
}

proc _enumStatesCountsToOffsets(counts, out _totalCounts) {
  const numChunks = counts.shape[0];

  var offsets : [0 ..# numChunks, 0 ..# numLocales] int;
  var totalCounts : [0 ..# numLocales] int;
  foreach localeIdx in 0 ..# numLocales {
    var total = 0;
    for chunkIdx in 0 ..# numChunks {
      offsets[chunkIdx, localeIdx] = total;
      total += counts[chunkIdx, localeIdx];
    }
    totalCounts[localeIdx] = total;
  }

  _totalCounts = totalCounts;
  return offsets;
}

proc _enumStatesMakeMasksOffsets(counts) {
  const numChunks = counts.shape[0];
  var masksOffsets : [0 ..# numChunks + 1] int;
  var total = 0;
  for chunkIdx in 0 ..# numChunks {
    masksOffsets[chunkIdx] = total;
    total += (+ reduce counts[chunkIdx, ..]);
  }
  masksOffsets[numChunks] = total;
  return masksOffsets;
}

proc _enumStatesMakeMasks(counts, totalCounts) {
  const numChunks = counts.shape[0];
  const numMasks = + reduce totalCounts;
  const masksBox = {0 ..# numMasks};
  const masksDom = masksBox dmapped Block(boundingBox=masksBox);
  var masks : [masksDom] uint(8);
  return masks;
}

proc _enumStatesPrepareDescriptors(masks, counts, totalCounts) {
  const numChunks = counts.shape[0];
  const masksOffsets = _enumStatesMakeMasksOffsets(counts);
  const numMasks = masksOffsets[numChunks];
  const ref masksDom = masks.domain;

  var masksPtrs : [0 ..# numLocales] c_ptr(uint(8));
  for loc in Locales {
    ref x = masks[masks.localSubdomain(loc).low];
    masksPtrs[loc.id] = __primitive("_wide_get_addr", x):c_ptr(uint(8));
  }

  var masksDescriptors : [0 ..# numChunks] (int, c_ptr(uint(8)));
  for i in 0 ..# numChunks {
    const firstIndex = masksOffsets[i];
    const lastIndex = if i < numChunks - 1
                        then masksOffsets[i + 1] - 1
                        else numMasks - 1;
    const loc = masksDom.dist.dsiIndexToLocale(firstIndex);
    // The range firstIndex .. lastIndex crosses the locale boundary in masks
    // array.
    if loc != masksDom.dist.dsiIndexToLocale(lastIndex) then
      masksDescriptors[i] = (firstIndex, nil);
    else {
      const ptrOffset = firstIndex - masksDom.localSubdomain(loc).low;
      masksDescriptors[i] = (loc.id, masksPtrs[loc.id] + ptrOffset);
    }
  }
  return masksDescriptors;
}

proc _enumStatesDistribute(const ref buckets, ref masks,
                           const ref counts,
                           const ref offsets,
                           const ref basisStatesPtrs,
                           const ref masksDescriptors) {
  var copyTimes : [0 ..# numLocales] real;
  var maskCopyTimes : [0 ..# numLocales] real;

  var distributeTimer = new Timer();
  distributeTimer.start();

  const numChunks = counts.shape[0];
  coforall loc in Locales do on loc {
    const mySubdomain = buckets.localSubdomain();
    const myCounts : [0 ..# numChunks, 0 ..# numLocales] int = counts;
    const myOffsets : [0 ..# numChunks, 0 ..# numLocales] int = offsets;
    const myDestPtrs : [0 ..# numLocales] c_ptr(uint(64)) = basisStatesPtrs;
    const myMasksDescriptors : [0 ..# numChunks] (int, c_ptr(uint(8))) = masksDescriptors;

    var myCopyTime : atomic real;
    var myMaskCopyTime : atomic real;
    forall bucketIdx in mySubdomain {
      const ref (myStates, myMasks) = buckets.localAccess(bucketIdx);

      // Distributing states
      const copyTime =
        permuteBasedOnMasks(myStates.size,
                            c_const_ptrTo(myMasks[0]),
                            c_const_ptrTo(myStates[0]),
                            myCounts[bucketIdx, ..],
                            myOffsets[bucketIdx, ..],
                            myDestPtrs);
      myCopyTime.add(copyTime, memoryOrder.relaxed);

      // Distributing masks
      var timer = new Timer();
      timer.start();
      const (localeOrOffset, targetPtr) = myMasksDescriptors[bucketIdx];
      if targetPtr != nil {
        const targetLocaleIdx = localeOrOffset;
        PUT(c_const_ptrTo(myMasks._arr[0]), targetLocaleIdx, targetPtr,
            myMasks.size:c_size_t * c_sizeof(uint(8)));
      }
      else {
        const targetOffset = localeOrOffset;
        masks[targetOffset ..# myMasks.size] = myMasks.toArray();
      }
      timer.stop();
      myMaskCopyTime.add(timer.elapsed(), memoryOrder.relaxed);
    }

    copyTimes[loc.id] = myCopyTime.read();
    maskCopyTimes[loc.id] = myMaskCopyTime.read();
  }
  distributeTimer.stop();

  return (distributeTimer.elapsed(), copyTimes, maskCopyTimes);
}

proc enumerateStates(ranges : [] range(uint(64)), const ref basis : Basis, out _masks) {
  var timer = new Timer();
  timer.start();

  // We distribute ranges among locales using Cyclic distribution to ensure
  // an even workload. For each range, a vector of basis states and a vector of
  // masks is computed. Masks indicate on which locale a basis state should live.
  const numChunks = ranges.size;
  var buckets = _enumStatesMakeBuckets(numChunks);

  // How many states coming from a certain chunk live on a certain locale.
  // Each chunk computes its own row in parallel and then does a remote PUT here.
  var countsTimer = new Timer();
  countsTimer.start();
  const counts = _enumStatesComputeCounts(buckets, ranges, basis);
  countsTimer.stop();

  // Transform counts into offsets such that each task known where to copy data.
  // Total counts tell us how much memory we have to allocate on each locale.
  const totalCounts;
  const offsets = _enumStatesCountsToOffsets(counts, totalCounts);

  // Allocate space for states
  var basisStates = new BlockVector(uint(64), totalCounts, distribute=true);
  const basisStatesPtrs = basisStates._dataPtrs;
  // Allocate space for masks
  var masks = _enumStatesMakeMasks(counts, totalCounts);
  var masksDescriptors =
    _enumStatesPrepareDescriptors(masks, counts, totalCounts);

  // Distribute buckets to basisStates and maks
  const (distributeTime, copyTimes, maskCopyTimes) =
    _enumStatesDistribute(buckets, masks, counts, offsets, basisStatesPtrs, masksDescriptors);

  timer.stop();
  logDebug("enumerateStatesNew took ", timer.elapsed(), "\n",
           "  ├─ ", countsTimer.elapsed(), " were spent computing counts\n",
           "  └─ ", distributeTime, " shuffling stuff around\n",
           "      ├─ ", copyTimes, " copying states\n",
           "      └─ ", maskCopyTimes, " copying masks");
  _masks = masks;
  return basisStates;
}
proc enumerateStates(globalRange : range(uint(64)), const ref basis : Basis,
                     out masks, numChunks : int = enumerateStatesNumChunks) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  const ranges = determineEnumerationRanges(globalRange, numChunks, isHammingWeightFixed);
  return enumerateStates(ranges, basis, masks);
}
proc enumerateStates(const ref basis : Basis, out masks,
                     numChunks : int = enumerateStatesNumChunks) {
  const lower = basis.minStateEstimate();
  const upper = basis.maxStateEstimate();
  return enumerateStates(lower .. upper, basis, masks, numChunks);
}


export proc ls_chpl_enumerate_representatives(p : c_ptr(ls_hs_basis),
                                              lower : uint(64),
                                              upper : uint(64),
                                              dest : c_ptr(chpl_external_array)) {
  logDebug("ls_chpl_enumerate_representatives ...");
  const basis = new Basis(p, owning=false);
  // var rs = localEnumerateRepresentatives(basis, lower, upper);
  const masks;
  var rs = enumerateStates(basis, masks);
  // ref v = rs[here];
  var v = rs[here];
  // writeln(v.type:string);
  // writeln(getExternalArrayType(v):string);
  // writeln(v._value.isDefaultRectangular():string);
  dest.deref() = convertToExternalArray(v);
  logDebug("Done!");
}

proc determineMergeBounds(const ref basisStates : [] uint(64), numChunks : int) {
  assert(basisStates.size >= numChunks);
  assert(numChunks > 1);
  var bounds : [0 ..# numChunks] uint(64);
  for (chunk, bound) in zip(chunks(0 ..# basisStates.size, numChunks), bounds) {
    bound = basisStates[chunk.high];
  }
  bounds[bounds.domain.low] = 0;
  bounds[bounds.domain.high] = ~(0:uint(64));
  return bounds;
}
// inline proc determineMergeBounds(const ref basisStates : Vector(uint(64)), numChunks : int) {
//   return determineMergeBounds(basisStates.toArray(), numChunks);
// }

proc mergeBoundsToIndexRanges(const ref basisStates : [] uint(64),
                              const ref bounds : [] uint(64)) {
  assert(bounds[bounds.domain.low] <= basisStates[basisStates.domain.low]);
  assert(bounds[bounds.domain.high] >= basisStates[basisStates.domain.high]);
  var ranges : [bounds.domain] range(int);
  for (r, b, i) in zip(ranges, bounds, 0..) {
    var loIdx = if i == 0 then 0 else ranges[i - 1].high + 1;
    var (found, hiIdx) = binarySearch(basisStates, b, lo=loIdx);
    if !found then
      hiIdx -= 1;
    r = loIdx .. hiIdx;
    // if hiIdx - loIdx >= 0 {
    //   if i != 0 {
    //     assert(basisStates[loIdx] > bounds[i - 1]);
    //   }
    //   assert(basisStates[hiIdx] <= b);
    // }
  }
  assert(ranges[bounds.size - 1].high == basisStates.size - 1);
  return ranges;
}
// inline proc mergeBoundsToIndexRanges(const ref basisStates : Vector(uint(64)),
//                                      const ref bounds : [] uint(64)) {
//   return mergeBoundsToIndexRanges(basisStates.toArray(), bounds);
// }

proc determineMergeRanges(const ref basisStates : BlockVector(uint(64), 1), numChunks : int) {
  // logDebug("determineMergeRanges: numStates=" + (+ reduce basisStates._counts):string);
  const bounds = determineMergeBounds(basisStates[here], numChunks);
  var ranges : [basisStates.outerDom] [0 ..# bounds.size] range(int);
  forall (r, i) in zip(ranges, 0 ..) {
    const localBounds = bounds;
    r = mergeBoundsToIndexRanges(basisStates.getBlock(i), localBounds);
  }
  return ranges;
}

proc mergeRangesToOffsets(const ref ranges) {
  const numRanges = ranges[0].size;
  var offsets : [0 ..# numRanges + 1] int;
  offsets[0] = 0;
  for i in 0 ..# numRanges {
    offsets[i + 1] = offsets[i] + (+ reduce [k in ranges.domain] ranges[k][i].size);
  }
  return offsets;
}

iter parallelMergeRangesIterator(basisStates : BlockVector(uint(64), 1), numChunks : int) {
  halt("not implemented, but required for compilation :(");
}

iter parallelMergeRangesIterator(param tag: iterKind, basisStates : BlockVector(uint(64), 1),
                                 in numChunks : int)
    where tag == iterKind.standalone {
  const maxNumChunks = min reduce basisStates.counts;
  if numChunks > maxNumChunks then
    numChunks = maxNumChunks;

  var timer = new Timer();
  timer.start();
  const ranges = determineMergeRanges(basisStates, numChunks);
  const offsets = mergeRangesToOffsets(ranges);
  const numRanges = offsets.size - 1;
  const numStates = offsets[numRanges];
  // logDebug("numStates=" + numStates:string);

  const space = {0 ..# numStates};
  const dom = space dmapped Block(boundingBox=space);
  timer.stop();
  logDebug("parallelMergeRangesIterator spent ", timer.elapsed(), " setting up");

  coforall rangeIdx in 0 ..# numRanges {
    const loc = dom.dist.dsiIndexToLocale(offsets[rangeIdx]);
    on loc {
      const localRange = offsets[rangeIdx] .. offsets[rangeIdx + 1] - 1;
      const localChunks = [i in LocaleSpace] ranges[i][rangeIdx];
      // logDebug("yielding: (" + localRange:string + ", " + localChunks:string + ")");
      yield (localRange, localChunks);
    }
  }
}

proc _hashedToBlockComputeCounts(const ref masks : [] ?i, numChunks : int) {
  var counts : [0 ..# numLocales, 0 ..# numChunks, 0 ..# numLocales] int;
  const countsPtr = c_ptrTo(counts[counts.domain.low]);
  coforall loc in Locales do on loc {
    const mySubdomain = masks.localSubdomain(loc);
    const myRanges : [0 ..# numChunks] range(int) =
      chunks(mySubdomain.dim(0), numChunks);
    // logDebug(myRanges, " and numChunks=", numChunks);
    // assert(myRanges.size == numChunks);
    var myCounts : [0 ..# numChunks, 0 ..# numLocales] int;
    forall (r, chunkIdx) in zip(myRanges, 0 ..# numChunks) {
      foreach key in masks.localAccess(r) {
        myCounts[chunkIdx, key:int] += 1;
      }
    }
    const myCountsPtr = c_const_ptrTo(myCounts[myCounts.domain.low]);
    const destPtr = countsPtr + loc.id * numChunks * numLocales;
    const destSize = (numChunks * numLocales):c_size_t * c_sizeof(int);
    PUT(myCountsPtr, 0, destPtr, destSize);
  }
  return counts;
}

proc _hashedToBlockComputeSrcOffsets(counts) {
  const numChunks = counts.shape[1];
  var offsets : [0 ..# numLocales, 0 ..# numChunks, 0 ..# numLocales] int;
  forall localeIdx in 0 ..# numLocales {
    var total = 0;
    for i in 0 ..# numLocales {
      for j in 0 ..# numChunks {
        offsets[i, j, localeIdx] = total;
        total += counts[i, j, localeIdx];
      }
    }
  }
  return offsets;
}

proc _hashedToBlockMakeDestArray(arr, masks) {
  param rank = arr.innerRank;
  const batchSize = if rank == 1 then 1
                                 else arr.innerDom.shape[0];
  const boundingBox = if rank == 1 then {0 ..# masks.size}
                                   else {0 ..# batchSize, 0 ..# masks.size};
  const targetLocales = if rank == 1 then Locales
                                     else reshape(Locales, {0 ..# 1, 0 ..# numLocales});
  const destArrDom = boundingBox dmapped Block(boundingBox, targetLocales);
  var destArr : [destArrDom] arr.eltType;
  return destArr;
}

proc _hashedToBlockNumChunks(masks, numChunks) {
  const minChunkSize = min reduce [loc in Locales] masks.localSubdomain(loc).size;
  return min(numChunks, minChunkSize);
}

proc arrFromHashedToBlock(const ref arr, const ref masks,
                          numChunks = _hashedToBlockNumChunks(masks,
                                        statesFromHashedToBlockNumChunks)) {
  var timer = new Timer();
  var countsTimer = new Timer();
  var distributeTimer = new Timer();

  timer.start();

  countsTimer.start();
  const counts = _hashedToBlockComputeCounts(masks, numChunks);
  countsTimer.stop();

  const srcOffsets = _hashedToBlockComputeSrcOffsets(counts);

  // const destArrBox = {0 ..# masks.size};
  // const destArrDom = destArrBox dmapped Block(boundingBox=destArrBox);
  var destArr = _hashedToBlockMakeDestArray(arr, masks);
  // : [destArrDom] arr.eltType;

  const countsPtr = c_const_ptrTo(counts[counts.domain.low]);
  const srcOffsetsPtr = c_const_ptrTo(srcOffsets[srcOffsets.domain.low]);
  const mainLocaleIdx = here.id;
  const arrPtrsPtr = c_const_ptrTo(arr._dataPtrs[0]);
  param rank = arr.innerRank;
  const batchSize = if rank == 1 then 1
                                 else arr.innerDom.shape[0];
  const batchStride = if rank == 1 then 0
                                   else arr.innerDom.shape[1];

  distributeTimer.start();
  coforall loc in Locales do on loc {
    const myDestSubdomain = destArr.localSubdomain();
    const myMasksSubdomain = masks.localSubdomain();
    // const mySize = myMasksSubdomain.size;
    const myRanges : [0 ..# numChunks] range(int) =
      chunks(myMasksSubdomain.dim(0), numChunks);

    var myCounts : [0 ..# numChunks, 0 ..# numLocales] int;
    var mySrcOffsets : [myCounts.domain] int;
    var myArrPtrs : [0 ..# numLocales] c_ptr(arr.eltType);
    cobegin {
      GET(c_ptrTo(myCounts[myCounts.domain.low]), mainLocaleIdx,
          countsPtr + loc.id * (numChunks * numLocales),
          (numChunks * numLocales):c_size_t * c_sizeof(int));
      GET(c_ptrTo(mySrcOffsets[mySrcOffsets.domain.low]), mainLocaleIdx,
          srcOffsetsPtr + loc.id * (numChunks * numLocales),
          (numChunks * numLocales):c_size_t * c_sizeof(int));
      GET(c_ptrTo(myArrPtrs[myArrPtrs.domain.low]), mainLocaleIdx,
          arrPtrsPtr, numLocales:c_size_t * c_sizeof(myArrPtrs.eltType));
    }
   
    for batchIdx in 0 ..# batchSize {
      forall (r, chunkIdx) in zip(myRanges, 0 ..# numChunks) {
        var maxCount = max reduce myCounts[chunkIdx, ..];
        var srcArr : [0 ..# numLocales, 0 ..# maxCount] arr.eltType;
        for srcLocaleIdx in 0 ..# numLocales {
          const n = myCounts[chunkIdx, srcLocaleIdx];
          const k = mySrcOffsets[chunkIdx, srcLocaleIdx];
          GET(c_ptrTo(srcArr[srcLocaleIdx, 0]), srcLocaleIdx,
              myArrPtrs[srcLocaleIdx] + batchIdx * batchStride + k,
              n:c_size_t * c_sizeof(srcArr.eltType));
        }

        var written : [0 ..# numLocales] int;
        for destIdx in r {
          const key = masks.localAccess(destIdx):int;
          ref offset = written[key];
          if rank == 1 then
            destArr.localAccess(destIdx) = srcArr[key, offset];
          else
            destArr.localAccess(batchIdx, destIdx) = srcArr[key, offset];
          offset += 1;
        }
      }
    }


  }
  distributeTimer.stop();

  timer.stop();
  logDebug("arrFromHashedToBlock took ", timer.elapsed(), "\n",
           "  ├─ ", countsTimer.elapsed(), " in computing counts\n",
           "  └─ ", distributeTimer.elapsed(), " in merging arrays");
  return destArr;
}


/*
proc statesFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1)) {
  const numStates = + reduce basisStates.counts; // + reduce basisStates._counts;
  const space = {0 ..# numStates};
  const dom = space dmapped Block(boundingBox=space);
  var blockBasisStates : [dom] uint(64);
  var blockBasisStatesPtrs : [0 ..# numLocales] c_ptr(uint(64));
  coforall loc in Locales do on loc {
    const mySubdomain = blockBasisStates.localSubdomain();
    blockBasisStatesPtrs[loc.id] = c_ptrTo(blockBasisStates[mySubdomain.low]);
  }
  const blockBasisStatesPtrsPtr = c_ptrTo(blockBasisStatesPtrs);

  const basisStatesDataPtrsPtr = c_ptrTo(basisStates._dataPtrs[0]);

  forall (currentRange, currentChunks) in
      parallelMergeRangesIterator(basisStates, statesFromHashedToBlockNumChunks) {
    var timer = new Timer();

    // timer.start();
    const currentCounts = [chunk in currentChunks] chunk.size;
    var localBasisStates = new BlockVector(uint(64), currentCounts, distribute=false);
    // timer.stop();
    // logDebug("statesFromHashedToBlock spent ", timer.elapsed(), " setting up");

    // Get data from remote
    timer.clear();
    timer.start();
    var basisStatesDataPtrs : [0 ..# numLocales] c_ptr(uint(64));
    if kUseLowLevelComm then
      GET(c_ptrTo(basisStatesDataPtrs[0]), 0, basisStatesDataPtrsPtr,
          numLocales:c_size_t * c_sizeof(c_ptr(uint(64))));
    for (i, chunk) in zip(0 ..# numLocales, currentChunks) {
      if kUseLowLevelComm then
        if chunk.size > 0 {
          const destPtr = c_ptrTo(localBasisStates[i, 0]);
          const destSize = chunk.size:c_size_t * c_sizeof(uint(64));
          const srcPtr = basisStatesDataPtrs[i] + chunk.low;
          GET(destPtr, i, srcPtr, destSize);
        }
      else
        localBasisStates._locBlocks[i][0 ..# chunk.size] = basisStates._locBlocks[i][chunk];
    }
    timer.stop();
    logDebug("statesFromHashedToBlock spent ", timer.elapsed(), " downloading");

    // Perform local merge
    timer.clear();
    timer.start();
    var mergedStates : [0 ..# currentRange.size] uint(64) = noinit;
    local {
      for ((chunkIndex, indexInChunk), offset) in zip(kMergeIndices(localBasisStates), 0..) {
        mergedStates[offset] = localBasisStates[chunkIndex, indexInChunk];
      }
    }
    timer.stop();
    logDebug("statesFromHashedToBlock spent ", timer.elapsed(), " merging");
    
    // Copy to remote
    timer.clear();
    timer.start();
    var blockBasisStatesPtrs : [0 ..# numLocales] c_ptr(uint(64));
    if kUseLowLevelComm then
      GET(c_ptrTo(blockBasisStatesPtrs[0]), 0, basisStatesDataPtrsPtr,
          numLocales:c_size_t * c_sizeof(c_ptr(uint(64))));
    // const targetLocale =
    //   blockBasisStates.domain.dist.dsiIndexToLocale(currentRange.low);
    // const targetSubdomain = blockBasisStates.localSubdomain(targetLocale);
    // if targetSubdomain.contains(currentRange) {
    //   
    // }
    // else {
      blockBasisStates[currentRange] = mergedStates;
    // }
    timer.stop();
    logDebug("statesFromHashedToBlock spent ", timer.elapsed(), " uploading");
  }

  return blockBasisStates;
}

proc makeBlockDomainForVectors(numVectors, numStates) {
  const boundingBox = {0 ..# numVectors, 0 ..# numStates};
  const targetLocales = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
  const dom : domain(2) dmapped Block(boundingBox=boundingBox, targetLocales=targetLocales) =
    boundingBox;
  return dom;
}

proc vectorsFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1),
                              const ref vectors : BlockVector(?eltType, 2)) {
  const numVectors = vectors.innerDom.dim(0).size;
  const numStates = + reduce basisStates.counts;
  const dom = makeBlockDomainForVectors(numVectors, numStates);
  var blockVectors : [dom] eltType;

  forall (currentRange, currentChunks) in
      parallelMergeRangesIterator(basisStates, statesFromHashedToBlockNumChunks) {
    const currentCounts = [chunk in currentChunks] chunk.size;
    var localBasisStates = new BlockVector(uint(64), currentCounts);
    var localVectors = new BlockVector(eltType, numVectors, currentCounts);
    // Get data from remote
    forall (i, chunk) in zip(LocaleSpace, currentChunks) {
      localBasisStates._data[i][0 ..# chunk.size] = basisStates._data[i][chunk];
      localVectors._data[i][.., 0 ..# chunk.size] = vectors._data[i][.., chunk];
    }

    // Perform local merge
    var mergedVectors : [0 ..# numVectors, 0 ..# currentRange.size] eltType = noinit;
    for ((chunkIndex, indexInChunk), offset) in zip(kMergeIndices(localBasisStates), 0..) {
      foreach k in 0 ..# numVectors {
        mergedVectors[k, offset] = localVectors._data[chunkIndex][k, indexInChunk];
      }
    }
    
    // Copy to remote
    blockVectors[.., currentRange] = mergedVectors;
  }

  return blockVectors;
}
*/

/*
proc statesAndVectorsFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1),
                                       const ref vectors : BlockVector(?eltType, 2)) {
  const ranges = determineMergeRanges(basisStates, statesFromHashedToBlockNumChunks);
  const offsets = mergeRangesToOffsets(ranges);
  const numRanges = offsets.size - 1;
  const numStates = offsets[numRanges];
  const numVectors = vectors._innerDom.dim(0).size;

  const statesDomain = {0 ..# numStates} dmapped Block(LocaleSpace);
  const vectorsDomain = makeBlockDomainForVectors(numVectors, numStates);
  var blockBasisStates : [statesDomain] uint(64);
  var blockVectors : [vectorsDomain] eltType;

  coforall rangeIdx in 0 ..# numRanges {
    const loc = statesDomain.dist.dsiIndexToLocale(offsets[rangeIdx]);
    on loc {
      const localRanges = [i in LocaleSpace] ranges[i][rangeIdx];
      const localCounts = [r in localRanges] r.size;
      var localBasisStates = new BlockVector(uint(64), localCounts);
      var localVectors = new BlockVector(eltType, numVectors, localCounts);
      forall (i, r) in zip(LocaleSpace, localRanges) {
        // ref v = localBasisStates.getBlock(i);
        localBasisStates._data[i][0 ..# r.size] = basisStates._data[i][r];
        localVectors._data[i][.., 0 ..# r.size] = vectors._data[i][.., r];
      }

      var mergedSize = offsets[rangeIdx + 1] - offsets[rangeIdx];
      var mergedStates : [0 ..# mergedSize] uint(64) = noinit;
      var mergedVectors : [0 ..# numVectors, 0 ..# mergedSize] eltType = noinit;


      for ((chunkIndex, indexInChunk), offset) in zip(kMergeIndices(localBasisStates), 0..) {
        mergedStates[offset] = localBasisStates._data[chunkIndex][indexInChunk];
        foreach k in 0 ..# numVectors {
          mergedVectors[k, offset] = localVectors._data[chunkIndex][k, indexInChunk];
        }
      }

      blockBasisStates[offsets[rangeIdx] .. offsets[rangeIdx + 1] - 1] = mergedStates;
      blockVectors[.., offsets[rangeIdx] .. offsets[rangeIdx + 1] - 1] = mergedVectors;
    }
  }
  return (blockBasisStates, blockVectors);
}
*/

/*
private proc distributionCounts(const ref basisStates : [] uint(64)) {
  var histogram : [LocaleSpace] int;
  forall basisState in basisStates with (+ reduce histogram) {
    const localeIdx = localeIdxOf(basisState);
    histogram[localeIdx] += 1;
  }

  const dom = LocaleSpace dmapped Block(LocaleSpace);
  var histogramDist : [dom] int = histogram;
  return histogramDist;
}



proc statesFromBlockToHashed(const ref basisStates : [] uint(64)) {
  const counts = distributionCounts(basisStates);
  var hashedBasisStates = new BlockVector(uint(64), counts);

  coforall loc in Locales with (ref hashedBasisStates) do on loc {
    ref dest = hashedBasisStates[loc];
    var offset = 0;

    for r in chunks(0 ..# basisStates.size, statesFromBlockToHashedNumChunks) {
      // Download from remote
      var localBasisStates = basisStates[r];
      // Determine which states to copy
      var mask : [0 ..# r.size] bool = noinit;
      forall (x, shouldInclude) in zip(localBasisStates, mask) {
        shouldInclude = localeIdxOf(x) == loc.id;
      }
      // Perform the copy
      for (x, shouldInclude) in zip(localBasisStates, mask) {
        if shouldInclude {
          dest[offset] = x;
          offset += 1;
        }
      }
    }
    assert(offset == dest.size);
  }
  return hashedBasisStates;
}
*/

/*
proc countingSort(src : [] ?t, dest : [] t, extractKey) {
  var counts : c_array(int, 256);
  foreach x in src {
    counts[extractKey(x)] += 1;
  }

  var total = 0;
  for i in 0 ..# 256 {
    const oldCount = counts[i];
    counts[i] = total;
    total += oldCount;
  }

  for x in src {
    const key = extractKey(x);
    dest[counts[key]] = x;
    counts[key] += 1;
  }
}


proc sortByLocaleIdx(const ref basisStates : [] uint(64),
                     ref vectors : [] ?eltType)
    where basisStates.domain.rank == 1 &&
          vectors.domain.rank == 2 {

  const numStates = basisStates.size;
  const numVectors = vectors.dim(0).size;
  if vectors.dim(1).size != numStates then
    halt("incompatible dimensions in sortByLocaleIdx: " + numStates:string +
         " vs. " + vectors.dim(1).size:string);

  if numLocales > 256 then
    halt("We're currently using uint(8) to store locale index. " +
         "Running with more than 256 locales is not yet supported");
  var masks : [0 ..# numStates] uint(8) = noinit;
  foreach (mask, basisState) in zip(masks, basisStates) {
    mask = localeIdxOf(basisState):uint(8);
  }

  var counts : [0 ..# numLocales] int;
  foreach mask in masks {
    counts[mask:int] += 1;
  }

  var offsets : [0 ..# numLocales] int;
  var sum : int = 0;
  for k in 0 ..# numLocales {
    offsets[k] = sum;
    sum += counts[k];
  }

  for i in 0 ..# numStates {
    
  }
}
*/

/*
proc vectorsFromBlockToHashed(const ref basisStates : [] uint(64),
                              const ref vectors : [] ?eltType) {
  const counts = distributionCounts(basisStates);
  const numVectors = vectors.dim(0).size;
  var hashedVectors = new BlockVector(eltType, numVectors, counts);

  coforall loc in Locales with (ref hashedVectors) do on loc {
    ref dest = hashedVectors[loc];
    var offset = 0;

    for r in chunks(0 ..# basisStates.size, statesFromBlockToHashedNumChunks) {
      // Download from remote
      var localBasisStates = basisStates[r];
      var localVectors = vectors[.., r];
      // Determine which states to copy
      var mask : [0 ..# r.size] bool = noinit;
      forall (x, shouldInclude) in zip(localBasisStates, mask) {
        shouldInclude = localeIdxOf(x) == loc.id;
      }
      // Perform the copy
      for i in 0 ..# r.size {
        const shouldInclude = mask[i];
        if shouldInclude {
          foreach k in 0 ..# numVectors {
            dest[k, offset] = localVectors[k, r.low + i];
          }
          offset += 1;
        }
      }
    }
    assert(offset == dest.dim(1).size);
  }
  return hashedVectors;
}
*/


proc _blockToHashedLocaleCounts(const ref masks : [] ?i) where isIntegral(i) {
  var counts : [0 ..# numLocales, 0 ..# numLocales] int;
  const countsPtr = c_ptrTo(counts[counts.domain.low]);
  const mainLocaleIdx = here.id;
  coforall loc in Locales do on loc {
    const mySubdomain = masks.localSubdomain();
    // const ref myMasks = masks[masks.localSubdomain()];
    var myCounts : [0 ..# numLocales] int;
    forall key in masks.localAccess(mySubdomain) with (+ reduce myCounts) {
      myCounts[key:int] += 1;
    }

    const putOffset = loc.id * numLocales;
    const putSize = numLocales:c_size_t * c_sizeof(int);
    PUT(c_ptrTo(myCounts[0]), mainLocaleIdx, countsPtr + putOffset, putSize);
  }
  return counts;
}

proc _blockToHashedLocaleOffsets(counts) {
  var offsets : [0 ..# numLocales, 0 ..# numLocales] int;
  foreach destLocaleIdx in 0 ..# numLocales {
    var total = 0;
    for srcLocaleIdx in 0 ..# numLocales {
      const count = counts[srcLocaleIdx, destLocaleIdx];
      offsets[srcLocaleIdx, destLocaleIdx] = total;
      total += count;
    }
  }
  return offsets;
}

proc _blockToHashedTaskCounts(masksSize : int, masksPtr : c_ptr(?i),
                              numChunksPerLocale : int) {
  // if kVerboseComm then startVerboseCommHere();
  var counts : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  const ranges : [0 ..# numChunksPerLocale] range(int) =
    chunks(0 ..# masksSize, numChunksPerLocale);

  forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
    foreach i in r do
      counts[chunkIdx, masksPtr[i]:int] += 1;
  }
  return counts;
}

proc _blockToHashedTaskOffsets(counts, perLocaleOffsets) {
  const numChunksPerLocale = counts.shape[0];
  var offsets : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  for destLocaleIdx in 0 ..# numLocales {
    var total = perLocaleOffsets[destLocaleIdx];
    for chunkIdx in 0 ..# numChunksPerLocale {
        const count = counts[chunkIdx, destLocaleIdx];
        offsets[chunkIdx, destLocaleIdx] = total;
        total += count;
    }
  }
  return offsets;
}

/*
*/

/*
proc permuteSmall(arrSize : int, masks : c_ptr(int), arr : c_ptr(?eltType),
                  counts : [] int, destOffsets : [] int, destPtrs : [] c_ptr(eltType)) {
  if kVerboseComm then startVerboseCommHere();
  // var prepareTimer = new Timer();
  // var computeTimer = new Timer();
  var copyTimer = new Timer();

  // prepareTimer.start();
  var offsets : [0 ..# numLocales] int = prefixSum(counts);
  var src : [0 ..# arrSize] eltType = noinit;
  // prepareTimer.stop();

  // computeTimer.start();
  for i in 0 ..# arrSize {
    const key = masks[i];
    src[offsets[key]] = arr[i];
    offsets[key] += 1;
  }
  // computeTimer.stop();

  copyTimer.start();
  var i = 0;
  for localeIdx in 0 ..# numLocales {
    if counts[localeIdx] > 0 {
      const srcPtr = c_ptrTo(src[i]);
      const destPtr = destPtrs[localeIdx] + destOffsets[localeIdx];
      const destSize = counts[localeIdx]:c_size_t * c_sizeof(eltType);
      PUT(srcPtr, localeIdx, destPtr, destSize);
      i += counts[localeIdx];
    }
  }
  copyTimer.stop();

  if kVerboseComm then stopVerboseCommHere();
  return copyTimer.elapsed();
}
*/

proc _blockToHashedMakeDestArr(arr : [] ?eltType, counts)
    where arr.domain.rank == 1 {
  const destCounts = [i in 0 ..# numLocales] (+ reduce counts[.., i]);
  return new BlockVector(eltType, destCounts);
}
proc _blockToHashedMakeDestArr(arr : [] ?eltType, counts)
    where arr.domain.rank == 2 {
  const destCounts = [i in 0 ..# numLocales] (+ reduce counts[.., i]);
  return new BlockVector(eltType, arr.shape[0], destCounts);
}

proc _blockToHashedNumChunksPerLocale(masks, numChunks) {
  const minBlockSize = min reduce [loc in Locales] masks.localSubdomain(loc).size;
  const suggested = max(1, numChunks / numLocales);
  return min(minBlockSize, suggested);
}

proc arrFromBlockToHashed(const ref arr : [] ?eltType, const ref masks : [] ?i,
                          numChunksPerLocale : int =
                            _blockToHashedNumChunksPerLocale(masks,
                              statesFromBlockToHashedNumChunks)) {
  var timer = new Timer();
  var countsTimer = new Timer();
  var makeDestArrTimer = new Timer();
  var distributeTimer = new Timer();
  var permuteTime : [0 ..# numLocales] real;
  const permuteTimePtr = c_ptrTo(permuteTime[0]);
  timer.start();

  countsTimer.start();
  const perLocaleCounts = _blockToHashedLocaleCounts(masks);
  const perLocaleOffsets = _blockToHashedLocaleOffsets(perLocaleCounts);
  const perLocaleOffsetsPtr = c_const_ptrTo(perLocaleOffsets[perLocaleOffsets.domain.low]);
  countsTimer.stop();

  makeDestArrTimer.start();
  var destArr = _blockToHashedMakeDestArr(arr, perLocaleCounts);
  const destPtrsPtr = c_ptrTo(destArr._dataPtrs);
  makeDestArrTimer.stop();

  distributeTimer.start();
  param rank = arr.domain.rank;
  const batchSize = if rank == 1
                      then 1
                      else arr.shape[0];
  const batchStride = if rank == 1
                        then 0
                        else destArr.innerDom.shape[1];
  const mainLocaleIdx = here.id;
  coforall loc in Locales do on loc {
    // const _myPtr : c_ptr(eltType) = c_ptrTo(destArr[loc.id, destArr.innerDom.low]);
    // if kUseLowLevelComm then
    //   PUT(c_const_ptrTo(_myPtr), 0, destPtrsPtr + loc.id, c_sizeof(c_ptr(eltType)));
    // else
    //   destPtrs[loc.id] = _myPtr;

    // allLocalesBarrier.barrier();

    var myPerLocaleOffsets : [0 ..# numLocales] int;
    GET(c_ptrTo(myPerLocaleOffsets[0]), mainLocaleIdx,
        perLocaleOffsetsPtr + loc.id * numLocales,
        numLocales:c_size_t * c_sizeof(int));

    var myDestPtrs : [0 ..# numLocales] c_ptr(eltType);
    GET(c_ptrTo(myDestPtrs[0]), mainLocaleIdx, destPtrsPtr,
        numLocales:c_size_t * c_sizeof(c_ptr(eltType)));

    const mySubdomain = masks.localSubdomain();
    const myMasksPtr = c_const_ptrTo(masks.localAccess(mySubdomain.low));
    const myMasksSize = mySubdomain.size;

    const myPerTaskCounts =
      _blockToHashedTaskCounts(myMasksSize, myMasksPtr, numChunksPerLocale);
    const myPerTaskOffsets = _blockToHashedTaskOffsets(myPerTaskCounts, myPerLocaleOffsets);
    // const (perTaskCount, perTaskOffset) =
    //   getPerTaskCountAndOffset(myMasksSize, myMasksPtr, myPerLocaleOffset);
    const ranges : [0 ..# numChunksPerLocale] range(int) =
      chunks(0 ..# myMasksSize, numChunksPerLocale);

    // var computeTimer = new Timer();
    // computeTimer.start();
    // var permuteSmallPrepareTime : atomic real;
    // var permuteSmallComputeTime : atomic real;
    var myPermuteTime : atomic real;

    for batchIdx in 0 ..# batchSize {
      const myArrPtr = if rank == 1
                         then c_const_ptrTo(arr.localAccess(mySubdomain.low))
                         else c_const_ptrTo(arr.localAccess(batchIdx, mySubdomain.low));

      forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
        var myTimer = new Timer();
        myTimer.start();
        permuteBasedOnMasks(r.size,
                            myMasksPtr + r.low,
                            myArrPtr + r.low,
                            myPerTaskCounts[chunkIdx, ..],
                            myPerTaskOffsets[chunkIdx, ..],
                            myDestPtrs);
        myTimer.stop();
        myPermuteTime.add(myTimer.elapsed(), memoryOrder.relaxed);
        // proc permuteBasedOnMasks(arrSize : int, masks : c_ptr(?maskType), arr : c_ptr(?eltType),
        //                          counts : [] int, destOffsets : [] int,
        //                          destPtrs : [] c_ptr(eltType)) {
        // const _copy =
        //   permuteSmall(r.size,
        //                myMasksPtr + r.low,
        //                myArrPtr + r.low,
        //                perTaskCount[chunkIdx, ..],
        //                perTaskOffset[chunkIdx, ..],
        //                myDestPtrs);
        // permuteSmallPrepareTime.add(_prep);
        // permuteSmallComputeTime.add(_compute);
        // permuteSmallCopyTime.add(_copy, memoryOrder.relaxed);
      }
      myDestPtrs += batchStride;
    }

    const _myPermuteTime = myPermuteTime.read();
    PUT(c_const_ptrTo(_myPermuteTime), mainLocaleIdx,
        permuteTimePtr + loc.id,
        c_sizeof(real));
    // computeTimer.stop();
    // logDebug("arrFromBlockToHashed main loop took ", computeTimer.elapsed(), "\n",
    //          "                      from which ", permuteSmallCopyTime.read(),
    //          " were spent in remote PUTs (divide by the number of tasks because ",
    //          "they ran in parallel)");
  }
  distributeTimer.stop();

  timer.stop();
  logDebug("arrFromBlockToHashed took ", timer.elapsed(), "\n",
           "  ├─ ", countsTimer.elapsed(), " were spent computing counts\n",
           "  ├─ ", makeDestArrTimer.elapsed(), " allocating destArr\n",
           "  └─ ", distributeTimer.elapsed(), " in the main loop\n",
           "      └─ ", permuteTime, " in permuteBasedOnMasks");
  // if kVerboseComm then stopVerboseComm();
  return destArr;
}


}
