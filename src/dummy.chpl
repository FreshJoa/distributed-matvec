use CTypes;
use BlockDist;
use RangeChunk;
use CommDiagnostics;
use AllLocalesBarriers;

use Vector;
use FFI;
/*
record PartitionInfo {
  var _countOrOffset : int;
  var nextOffset : int;

  inline proc count ref { return _countOrOffset; }
  inline proc offset ref { return _countOrOffset; }
}

proc _partition(in first : c_ptr(?t), last : c_ptr(t), predicate) {
    while first != last {
      if !predicate(first) then
        break;
      first += 1;
    }
    if first == last then
      return first;

    var i = first + 1;
    while i != last {
      if predicate(i) {
        i.deref() <=> first.deref();
        first += 1;
      }
      i += 1;
    }
    return first;
}

proc skaByteSort(ref arr : [] ?eltType, extractKey) {
  var partitions : c_array(PartitionInfo, 256);
  foreach x in arr {
    partitions[extractKey(x)].count += 1;
  }

  var remainingPartitions : c_array(uint(8), 256);
  var numPartitions : int;
  var total : int;

  for i in 0 ..# 256 {
    const count = partitions[i].count;
    if count > 0 {
      partitions[i].offset = total;
      total += count;
      remainingPartitions[numPartitions] = i:uint(8);
      numPartitions += 1;
    }
    partitions[i].nextOffset = total;
  }

  var lastRemaining = remainingPartitions:c_ptr(uint(8)) + numPartitions;
  var endPartition = remainingPartitions:c_ptr(uint(8)) + 1;
  while lastRemaining - endPartition > 0 {
    record Predicate {
      inline proc this(p : c_ptr(uint(8))) : bool {
        ref beginOffset = partitions[p.deref()].offset;
        ref endOffset = partitions[p.deref()].nextOffset;
        if beginOffset == endOffset then
          return false;

        for i in beginOffset .. endOffset - 1 {
          const key = extractKey(arr[i]);
          const offset = partitions[key].offset;
          partitions[key].offset += 1;
          arr[i] <=> arr[offset];
        }
        return beginOffset != endOffset;
      }
    }
    lastRemaining = _partition(remainingPartitions:c_ptr(uint(8)),
                               lastRemaining,
                               new Predicate());

  }
}
*/

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

inline proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

inline proc localeIdxOf(basisState : uint(64)) : int {
  return (hash64_01(basisState) % numLocales:uint):int;
}

config const numChunksPerLocale = 3;
config const kUseLowLevelComm : bool = true;
config const kVerboseComm : bool = false;

inline proc GET(addr, node, rAddr, size) {
  __primitive("chpl_comm_get", addr, node, rAddr, size);
}

inline proc PUT(addr, node, rAddr, size) {
  __primitive("chpl_comm_put", addr, node, rAddr, size);
}

proc getPerLocaleCountAndOffset(masks : [] int) {
  if kVerboseComm then startVerboseCommHere();
  var perLocaleCount : [0 ..# numLocales, 0 ..# numLocales] int;
  const perLocaleCountPtr = c_ptrTo(perLocaleCount[0, 0]);
  coforall loc in Locales do on loc {
    const ref myMasks = masks[masks.localSubdomain()];
    var myCount : [0 ..# numLocales] chpl__processorAtomicType(int);
    forall key in myMasks {
      myCount[key].add(1, memoryOrder.relaxed);
    }

    if kUseLowLevelComm {
      const putOffset = loc.id * numLocales;
      const putSize = numLocales:c_size_t * c_sizeof(int);
      PUT(c_ptrTo(myCount[0]):c_ptr(int), 0, perLocaleCountPtr + putOffset, putSize);
    }
    else {
      perLocaleCount[loc.id, ..] = [i in 0 ..# numLocales] myCount[i].read();
    }
  }

  var perLocaleOffset : [0 ..# numLocales, 0 ..# numLocales] int;
  foreach destLocaleIdx in 0 ..# numLocales {
    var total = 0;
    for srcLocaleIdx in 0 ..# numLocales {
      const count = perLocaleCount[srcLocaleIdx, destLocaleIdx];
      perLocaleOffset[srcLocaleIdx, destLocaleIdx] = total;
      total += count;
    }
  }

  if kVerboseComm then stopVerboseCommHere();
  return (perLocaleCount, perLocaleOffset);
}

proc getPerTaskCountAndOffset(masksSize : int, masksPtr : c_ptr(int),
                              perLocaleOffset : [] int) {
  if kVerboseComm then startVerboseCommHere();
  var perTaskCount : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  var perTaskOffset : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  const ranges : [0 ..# numChunksPerLocale] range(int) =
    chunks(0 ..# masksSize, numChunksPerLocale);

  forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
    foreach i in r do
      perTaskCount[chunkIdx, masksPtr[i]] += 1;
  }

  for destLocaleIdx in 0 ..# numLocales {
    var total = perLocaleOffset[destLocaleIdx];
    for chunkIdx in 0 ..# numChunksPerLocale {
        const count = perTaskCount[chunkIdx, destLocaleIdx];
        perTaskOffset[chunkIdx, destLocaleIdx] = total;
        total += count;
    }
  }
  if kVerboseComm then stopVerboseCommHere();
  return (perTaskCount, perTaskOffset);
}

proc prefixSum(arr : [] ?eltType) {
  var sums : [0 ..# arr.size] eltType;
  var total : eltType = 0;
  for i in sums.domain {
    sums[i] = total;
    total += arr[i];
  }
  return sums;
}

proc permuteSmall(arrSize : int, masks : c_ptr(int), arr : c_ptr(?eltType),
                  counts : [] int, destOffsets : [] int, destPtrs : [] c_ptr(eltType)) {
  if kVerboseComm then startVerboseCommHere();
  var offsets : [0 ..# numLocales] int = prefixSum(counts);
  var src : [0 ..# arrSize] eltType;
  for i in 0 ..# arrSize {
    const key = masks[i];
    src[offsets[key]] = arr[i];
    offsets[key] += 1;
  }

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
  if kVerboseComm then stopVerboseCommHere();
}

proc _makeDestArr(arr : [] ?eltType, perLocaleCount)
    where arr.domain.rank == 1 {
  const dom = LocaleSpace dmapped Block(LocaleSpace, Locales);
  const destCounts = [i in 0 ..# numLocales] (+ reduce perLocaleCount[.., i]);
  return new BlockVector(eltType, destCounts, dom);
}
proc _makeDestArr(arr : [] ?eltType, perLocaleCount)
    where arr.domain.rank == 2 {
  const dom = LocaleSpace dmapped Block(LocaleSpace, Locales);
  const destCounts = [i in 0 ..# numLocales] (+ reduce perLocaleCount[.., i]);
  return new BlockVector(eltType, arr.shape[0], destCounts, dom);
}

proc arrFromBlockToHashed(masks : [] int, arr : [] ?eltType) {
  const (perLocaleCount, perLocaleOffset) = getPerLocaleCountAndOffset(masks);
  const perLocaleOffsetPtr = c_const_ptrTo(perLocaleOffset[0, 0]);
  var destArr = _makeDestArr(arr, perLocaleCount);
  var destPtrs : [0 ..# numLocales] c_ptr(eltType);
  const destPtrsPtr = c_const_ptrTo(destPtrs);

  if kVerboseComm then startVerboseComm();
  const batchSize = if arr.domain.rank == 1 then 1 else arr.shape[0];
  const batchIncrement = if arr.domain.rank == 1
                           then 0 else destArr._innerDom.shape[1];
  coforall loc in Locales do on loc {
    const _myPtr : c_ptr(eltType) = c_ptrTo(destArr._data[loc.id][destArr._innerDom.low]);
    if kUseLowLevelComm then
      PUT(c_const_ptrTo(_myPtr), 0, destPtrsPtr + loc.id, c_sizeof(c_ptr(eltType)));
    else
      destPtrs[loc.id] = _myPtr;

    allLocalesBarrier.barrier();

    var myPerLocaleOffset : [0 ..# numLocales] int = noinit;
    if kUseLowLevelComm then
      GET(c_ptrTo(myPerLocaleOffset[0]), 0, perLocaleOffsetPtr + loc.id * numLocales,
          numLocales:c_size_t * c_sizeof(int));
    else
      myPerLocaleOffset = perLocaleOffset[loc.id, ..];

    var myDestPtrs : [0 ..# numLocales] c_ptr(eltType) = noinit;
    if kUseLowLevelComm then
      GET(c_ptrTo(myDestPtrs[0]), 0, destPtrsPtr,
          numLocales:c_size_t * c_sizeof(c_ptr(eltType)));
    else
      myDestPtrs = destPtrs;

    const mySubdomain = masks.localSubdomain();
    const myMasksPtr = c_const_ptrTo(masks[mySubdomain.low]);
    const myMasksSize = mySubdomain.size;

    const (perTaskCount, perTaskOffset) =
      getPerTaskCountAndOffset(myMasksSize, myMasksPtr, myPerLocaleOffset);
    const ranges : [0 ..# numChunksPerLocale] range(int) =
      chunks(0 ..# myMasksSize, numChunksPerLocale);

    for batchIdx in 0 ..# batchSize {
      const myArrPtr = if arr.domain.rank == 1
                         then c_const_ptrTo(arr[mySubdomain.low])
                         else c_const_ptrTo(arr[batchIdx, mySubdomain.low]);

      forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
        permuteSmall(r.size,
                     myMasksPtr + r.low,
                     myArrPtr + r.low,
                     perTaskCount[chunkIdx, ..],
                     perTaskOffset[chunkIdx, ..],
                     myDestPtrs);
      }
      myDestPtrs += batchIncrement;
    }
  }
  if kVerboseComm then stopVerboseComm();
  return destArr;
}

/*
proc generateOffsetsFromLocaleIndices(masks : [] int) {
  startVerboseComm();
  // masks is a block-distributed array of locale indices
  var minorOffsets : [0 ..# numLocales, 0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  var majorCounts : [0 ..# numLocales, 0 ..# numLocales] int;
  const minorOffsetsPtr = c_ptrTo(minorOffsets[0, 0, 0]);
  const majorCountsPtr = c_ptrTo(majorCounts[0, 0]);

  coforall loc in Locales do on loc {
    var myMinorCounts : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
    var myMajorCounts : [0 ..# numLocales] int;
    const ref myMasks = masks[masks.localSubdomain()];
    local {
      const ranges : [0 ..# numChunksPerLocale] range(int) =
        chunks(myMasks.domain.low .. myMasks.domain.high, numChunksPerLocale);

      forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
        foreach key in myMasks.localAccess(r) {
          myMinorCounts[chunkIdx, key] += 1;
        }
      }

      foreach localeIdx in 0 ..# numLocales do
        myMajorCounts[localeIdx] = + reduce myMinorCounts[.., localeIdx];

      // Transform myMinorCounts into offsets
      for chunkIdx in 0 ..# numChunksPerLocale {
        var total = 0;
        for localeIdx in 0 ..# numLocales {
          const oldOffset = myMinorCounts[chunkIdx, localeIdx];
          myMinorCounts[chunkIdx, localeIdx] = total;
          total += oldOffset;
        }
      }

      // for localeIdx in 0 ..# numLocales {
      //   var total = 0;
      //   for chunkIdx in 0 ..# numChunksPerLocale {
      //     const oldOffset = localOffsets[chunkIdx, localeIdx];
      //     localOffsets[chunkIdx, localeIdx] = total;
      //     total += oldOffset;
      //   }
      //   localOffsets[numChunksPerLocale, localeIdx] = total;
      // }
    } // end local

    // We want to do
    //   globalOffsets[loc.id, .., ..] = localOffsets;
    // but this is currently inefficient
    {
      const destOffset = loc.id * numChunksPerLocale * numLocales;
      const destSize = numChunksPerLocale * numLocales;
      PUT(c_ptrTo(myMinorCounts[0, 0]),
          0, minorOffsetsPtr + destOffset,
          destSize:c_size_t * c_sizeof(int));
    }
    {
      const destOffset = loc.id * numLocales;
      const destSize = numLocales;
      PUT(c_ptrTo(myMajorCounts[0]),
          0, majorCountsPtr + destOffset,
          destSize:c_size_t * c_sizeof(int));
    }

    stopVerboseCommHere();

    // for chunkIdx in 0 ..# numChunksPerLocale {
    //   writeln(loc, ": ", ranges[chunkIdx], ": ", localOffsets[chunkIdx, ..]);
    // }
    // writeln(loc, ": _____: ", localOffsets[numChunksPerLocale, ..]);
  }

  // var majorOffsets : [0 ..# numLocales, 0 ..# numLocales] int;
  // for localeIdx in 0 ..# numLocales {
  //   var total = 0;
  //   for k in 0 ..# numLocales {
  //     const minorOffset = minorOffsets[k, numChunksPerLocale, localeIdx];
  //     majorOffsets[k, localeIdx] += total;
  //     total += minorOffset;
  //   }
  // }

  for loc in Locales {
    for chunkIdx in 0 ..# numChunksPerLocale {
      writeln(loc, ": ", minorOffsets[loc.id, chunkIdx, ..]);
    }
    // writeln(loc, ": ", minorOffsets[loc.id, numChunksPerLocale, ..]);
    writeln(loc, ": major counts: ", majorCounts[loc.id, ..]);
  }
  stopVerboseComm();
  return (majorCounts, minorOffsets);
}

proc permuteStatesBasedOnLocaleIndices(masks : [] int,
                                       basisStates : [] ?t) {
  const (majorCounts, minorOffsets) = generateOffsetsFromLocaleIndices(masks);
  var targetCounts = [i in 0 ..# numLocales] (+ reduce majorCounts[.., i]);
  writeln(targetCounts);

  const dom = LocaleSpace dmapped Block(LocaleSpace, Locales);
  var hashedStates = new BlockVector(int, targetCounts, dom);

  for loc in Locales do on loc {
    // const myMajorOffsets = majorOffsets[loc.id, ..];
    const myMinorOffsets = minorOffsets[loc.id, .., ..];
    const mySubdomain = masks.localSubdomain();
    const ref myMasks = masks[mySubdomain];
    const ref myStates = basisStates[mySubdomain];

    const ranges : [0 ..# numChunksPerLocale] range(int) =
      chunks(mySubdomain.low .. mySubdomain.high, numChunksPerLocale);

    for (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
      var localOffsets = myMinorOffsets[chunkIdx, ..];
      writeln(localOffsets);

      // for x in myMasks[r] {
      //   const key = localeIdxOf(x:uint);
      //   localOffsets[chunkIdx, key] += 1;
      // }
  // for x in src {
  //   const key = extractKey(x);
  //   dest[counts[key]] = x;
  //   counts[key] += 1;
  // }
    }

  }
  /*


  // startVerboseComm();
  // startCommDiagnostics();
  // stopCommDiagnostics();
  // for (info, loc) in zip(getCommDiagnostics(), Locales) {
  //   writeln(loc, ": ", info);
  // }
  // stopVerboseComm();


  */
}
*/


proc main() {
  const box = {0 ..# 20};
  const dom = box dmapped Block(box, Locales);
  var arr : [dom] uint = 0:uint ..# 20:uint;
  const masks = [x in arr] localeIdxOf(x);

  writeln(arr);
  writeln(masks);

  const box2 = {0 ..# 2, 0 ..# 20};
  const targetLocales2 = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
  const dom2 = box2 dmapped Block(box2, targetLocales2);
  var arr2 : [dom2] real;
  for i in 0 ..# 2 {
    arr2[i, ..] = (0 ..# 20):real;
  }
  arr2[1, ..] *= -1;

  // var counts : [0 ..# numLocales] int;
  // for key in masks {
  //   counts[key] += 1;
  // }

  const newArr = arrFromBlockToHashed(masks, arr2);
  for loc in Locales {
    writeln(loc, ": ", newArr[loc]);
  }


  // generateVariousOffsets(masks);

  // generateOffsetsFromLocaleIndices(arr);
  // permuteStatesBasedOnLocaleIndices(masks, arr);

  // record Predicate {
  //   inline proc this(i : c_ptr(int)) : bool {
  //     return i.deref() % 2 == 0;
  //   }
  // }

  // const first = c_ptrTo(arr[arr.domain.low]);
  // const last = first + arr.size;
  // const i = _partition(first, last, new Predicate());
  // writeln(i - first);

  // writeln(arr);
}
