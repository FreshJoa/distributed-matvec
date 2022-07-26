module Distribute {
  use BlockDist;
  use CommDiagnostics;
  // use BlockCycDist;
  // use CyclicDist;
  use states only hash64_01;
  use basis only BasisStates;
  use profiling;

  config const kDistributionMaskChunkSize = 10240;

  var _distributionMaskTime = new MeasurementTable("distributionMask");
  var _distributionMaskInnerTime = new MeasurementTable("distributionMask::inner copy");

  proc distributionMask(const ref states : [?D] uint(64)) {
    const n = numLocales;
    var mask : [D] int;
    forall (x, i) in zip(states, mask) {
      assert(x.locale == here, "x is on the wrong locale");
      assert(i.locale == here, "i is on the wrong locale");
      i = (hash64_01(x) % n:uint):int;
    }
    return mask;
  }

  proc distributionMask(const ref states : BasisStates,
                        chunkSize : int = kDistributionMaskChunkSize) {
    var __timer = getTimerFor(_distributionMaskTime);
    const totalCount = states.totalNumberStates();
    const boundingBox = {0 ..# totalCount};
    const D : domain(1) dmapped Block(boundingBox=boundingBox) = boundingBox;
    var mask : [D] int;
    forall (offset, indices) in kMergeIndicesChunked(states._representatives,
                                                     states._counts, chunkSize)
        with (ref _distributionMaskInnerTime) {
      var __timer1 = getTimerFor(_distributionMaskInnerTime);
      // if indices.locale != here { error("indices are on the wrong locale"); }
      // if offset.locale == here { error("offset is on the wrong locale"); }

      const localMask : [0 ..# indices.size] int = [i in indices.domain] indices[i][0];
      mask[offset ..# indices.size] = localMask;
    }
    // writeln(D);
    // for i in LocaleSpace {
    //   writeln(mask[mask.localSubdomain(Locales[i])]);
    // }
    return mask;
  }

  private proc distributionCounts(const ref mask : [] int) {
    var histogram : [LocaleSpace] int;
    forall i in mask with (+ reduce histogram) {
      histogram[i] += 1;
    }
    return histogram;
  }

  private proc maxInnerCount(const ref mask : [] int)
      where mask.hasSingleLocalSubdomain() {
    var globalCounts : [LocaleSpace] [LocaleSpace] int;
    coforall loc in Locales do on loc {
      const localIndices = mask.localSubdomain();
      const localCounts = distributionCounts(mask[localIndices]);
      globalCounts[loc.id] = localCounts;
    }
    return max reduce ([i in LocaleSpace] max reduce globalCounts[i]);
  }

  private inline proc innerDomain(const ref xs : [?D] ?eltType, n : int) : domain(1)
      where D.rank == 1 {
    return {0 ..# n};
  }
  private inline proc innerDomain(const ref xs : [?D] ?eltType, n : int) : domain(2)
      where D.rank == 2 {
    return {0 ..# D.dim(0).size, 0 ..# n};
  }

  var _distributeArrayTime = new MeasurementTable("distributeArray");
  var _distributeArrayDistributeTime = new MeasurementTable("distributeArray::distribute");
  var _distributeArrayRemoteCopiesTime = new MeasurementTable("distributeArray::remote copy");

  proc distributeArray(const ref array : [?OuterDomain] ?eltType,
                       const ref mask  : [?OuterDomain2] int)
      where isSubtype(OuterDomain.dist.type, Block) {
    assert(array.size > 0);
    var __timer = getTimerFor(_distributeArrayTime);
    param rank = array.domain.rank;
    const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);

    var __timer1 = getTimerFor(_distributeArrayDistributeTime);
    const TempInnerDomain = innerDomain(array, maxInnerCount(mask));
    var globalTemp : [OnePerLocale] [LocaleSpace] [TempInnerDomain] eltType;
    var globalTempCounts : [OnePerLocale] [LocaleSpace] int;
    // startVerboseComm();
    coforall loc in Locales do on loc {
      ref temp = globalTemp[loc.id];
      ref offsets = globalTempCounts[loc.id];
      // forall targetLocale in LocaleSpace {
        for i in mask.localSubdomain() do {
          const targetLocale = mask[i];
          // if mask[i] == targetLocale {
            if rank == 1 {
              temp[targetLocale][offsets[targetLocale]] = array.localAccess[i];
            }
            else if rank == 2 {
              foreach k in TempInnerDomain.dim(0) {
                temp[targetLocale][k, offsets[targetLocale]] = array.localAccess[k, i];
              }
            }
            else { compilerError("Oops"); }
            offsets[targetLocale] += 1;
          // }
        }
      // }
    }
    // stopVerboseComm();
    __timer1.stop();

    var __timer2 = getTimerFor(_distributeArrayRemoteCopiesTime);
    const globalFinalCounts : [OnePerLocale] int =
      [j in LocaleSpace] (+ reduce [i in OnePerLocale] globalTempCounts[i][j]);
    const maxCountPerLocale = max reduce globalFinalCounts;
    const FinalInnerDomain = innerDomain(array, maxCountPerLocale);
    var globalFinal : [OnePerLocale] [FinalInnerDomain] eltType;
    coforall loc in Locales do on loc {
      var offset : int = 0;
      for i in LocaleSpace {
        const n = globalTempCounts[i][loc.id];
        if rank == 1 {
          globalFinal[loc.id][offset ..# n] = globalTemp[i][loc.id][0 ..# n];
        }
        else if rank == 2 {
          globalFinal[loc.id][.., offset ..# n] = globalTemp[i][loc.id][.., 0 ..# n];
        }
        offset += n;
      }
    }
    __timer2.stop();

    return (globalFinal, globalFinalCounts);
  }

  inline proc distributeStates(const ref states) {
    return distributeStates(states, distributionMask(states));
  }
  proc distributeStates(const ref states : [?D] uint(64),
                        const ref mask : [?D2] int)
      where D.rank == 1 {
    var (chunks, counts) = distributeArray(states, mask); 
    return new BasisStates(chunks[0].size, chunks, counts);
  }

  proc distributeVectors(const ref vectors : [?D] ?eltType,
                         const ref mask : [?D2] int)
      where D.rank == 2 {
    var (chunks, counts) = distributeArray(vectors, mask); 
    return chunks;
  }
} // module Distribute
