module MatVec {
  module IO {
    use wrapper;
    use basis only OnePerLocale;
    use Search;
    use Time;
    use BlockDist;
    use Merge;
    use profiling;

    config const kMergeChunkSize : int = 5;
    // config const kLoadChunkSize : int = 10 * 1024 / numLocales;

    proc gatherImpl(param tag : int, const ref arrays, const ref indices : [?D] (int, int))
        where tag == 1 {
      type eltType = arrays[0].eltType;
      var combined : [0 ..# indices.size] eltType;
      forall ((arrayIndex, indexInArray), offset) in zip(indices, 0..) {
        combined[offset] = arrays[arrayIndex][indexInArray];
      }
      return combined;
    }
    proc gatherImpl(param tag : int, const ref arrays, const ref indices : [?D] (int, int))
        where tag == 2 {
      type eltType = arrays[0].eltType;
      const numberVectors = arrays[0].dim(0).size;
      var combined : [0 ..# numberVectors, 0 ..# indices.size] eltType;
      forall ((arrayIndex, indexInArray), offset) in zip(indices, 0..) {
        forall k in 0 ..# numberVectors {
          combined[k, offset] = arrays[arrayIndex][k, indexInArray];
        }
      }
      return combined;
    }

    inline proc gather(const ref arrays, const ref indices : [?D] (int, int))
        where D.rank == 1 {
      return gatherImpl(arrays[0].domain.rank, arrays, indices);
    }

    /* Merge states and write them to HDF5 file.
     *
     * :arg states: distributed array of sorted arrays
     * :type states: [] [] uint(64)
     *
     * :arg counts: specifies number of elements in each array of `states`
     * :type counts: [] int
     *
     * :arg filename: path to output HDF5 file.
     * :arg dataset:  path to dataset within HDF5 file.
     */
    proc mergeAndWriteStates(const ref states, const ref counts : [] int,
                             filename: string, dataset: string,
                             chunkSize : int = 5) {
      const totalCount = (+ reduce counts);
      createHDF5Dataset(filename, dataset, states[0].eltType, (totalCount,));
      forall (offset, indices) in kMergeIndicesChunked(states, counts, chunkSize) {
        var localStates = gather(states, indices);
        // writeln("[Chapel] writing ", localStates);
        writeHDF5Chunk(filename, dataset, (offset,), localStates);
      }
    }

    var _mergeAndWriteVectorsTime = new MeasurementTable("mergeAndWriteVectors");

    proc mergeAndWriteVectors(const ref states, const ref vectors, const ref counts : [] int,
                              filename: string, dataset: string,
                              chunkSize : int = 5) {
      var __timer = getTimerFor(_mergeAndWriteVectorsTime);
      const totalCount = (+ reduce counts);
      const numberVectors = vectors[0].dim(0).size;
      createHDF5Dataset(filename, dataset, vectors[0].eltType, (numberVectors, totalCount));
      forall (offset, indices) in kMergeIndicesChunked(states, counts, chunkSize) {
        var localVectors = gather(vectors, indices);
        // writeln("[Chapel] writing ", localStates);
        writeHDF5Chunk(filename, dataset, (0, offset), localVectors);
      }
    }

    var _loadStatesTime = new MeasurementTable("loadStates");

    proc loadStates(filename : string, dataset : string) {
      var __timer = getTimerFor(_loadStatesTime);
      const _shape = datasetShape(filename, dataset);
      assert(_shape.size == 1);
      const totalNumberStates = _shape[0];
     
      const D : domain(1) dmapped Block(LocaleSpace) = {0 ..# totalNumberStates};
      var states : [D] uint(64);
      coforall loc in Locales do on loc {
        const indices = states.localSubdomain();
        readHDF5Chunk(filename, dataset, (indices.low,), states[indices]);
      }
      return states;
    }

    var _loadVectorsTime = new MeasurementTable("loadVectors");

    proc loadVectors(filename : string, dataset : string, type eltType = real(64)) {
      var __timer = getTimerFor(_loadVectorsTime);
      const _shape = datasetShape(filename, dataset);
      assert(_shape.size == 2);
      const numberVectors = _shape[0];
      const totalNumberStates = _shape[1];
    
      const boundingBox = {0 ..# numberVectors, 0 ..# totalNumberStates};
      const targetLocales = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
      const D : domain(2) dmapped Block(boundingBox=boundingBox, targetLocales=targetLocales) =
        boundingBox;
      var vectors : [D] eltType;
      coforall loc in Locales do on loc {
        const indices = vectors.localSubdomain();
        readHDF5Chunk(filename, dataset, indices.low, vectors[indices]);
      }
      for loc in Locales {
        writeln("vectors[", loc.id, "]: ", vectors.localSubdomain(loc));
      }
      return vectors;
    }

    // proc mergeAndWriteVectors(const ref states, const ref vectors, const ref counts : [] int,
    //                           filename: string, dataset: string) {
    //   MergeAndWrite.mainLoop(counts, filename, (none, dataset), states, vectors);
    // }

    /*
    proc mergeAndWriteStates(const ref states, const ref counts : [] int,
                             filename: string, dataset: string) {
      // Prepare the output file
      const numberChunks = counts.size;
      const totalCount = (+ reduce counts);
      createHDF5Dataset(filename, dataset, uint(64), (totalCount,));

      // Split all states into smaller chunks
      const indexEdges = computeRanges(states, counts);
      const numberRanges = indexEdges.domain.dim(1).size - 1;
      const offsets = computeOffsets(indexEdges);
      // Compute sizes of those chunks
      // var offsets : [0 .. numberRanges] uint(64);
      // offsets[0] = 0;
      // for j in 1 .. numberRanges {
      //   const n = (+ reduce [i in indexEdges.dim(0)] indexEdges[i, j] - indexEdges[i, j - 1]);
      //   offsets[j] = offsets[j - 1] + n:uint;
      // }
      // Iterate over chunks
      // TODO: parallelize
      for k in 0 .. numberRanges - 1 {
        // Prepare for kmerge
        var localCounts = [i in indexEdges.dim(0)] indexEdges[i, k + 1] - indexEdges[i, k];
        // var maxLocalCount = max reduce localCounts;
        // Copy chunks from different locales to here
        (localStates, localArrays) = _createLocalCopyForMerge(states,
            localCounts, maxLocalCount, indexEdges[i, k], indexEdges[i, k + 1]);
        // var localStates : [0 .. states.dim(0).size - 1, 0 .. maxLocalCount - 1] uint(64);
        // for i in localStates.dim(0) {
        //   localStates[i, 0 .. localCounts[i] - 1] =
        //     states[i][indexEdges[i, k] .. indexEdges[i, k + 1] - 1];
        // }
        // Run kmerge locally

        // const localCount = offsets[k + 1] - offsets[k];
        // assert(localCount:int == + reduce localCounts);
        // TODO: not const because of unsupported const to non-const pointer casts ...
        var combinedLocalStates = _localMergeStates(localStates, localVectors, localCounts);

        // Save to file
        // TODO: make async
        assert(combinedLocalStates.size == (offsets[k + 1] - offsets[k]):int);
        writeHDF5Chunk(filename, dataset, (offsets[k],), combinedLocalStates);
        // var c_offset = offsets[k];
        // var c_shape = localCount;
        // ls_hs_hdf5_write_chunk_u64(filename.c_str(), dataset.c_str(),
        //   1, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(combinedLocalStates));
      }
    }

    proc mergeVectorsBy(const ref states, const ref counts : [] int, const ref vectors,
                        filename: string, dataset: string) {
      // Prepare the output file
      const totalCount = (+ reduce counts);
      assert(vectors[0].rank == 2);
      const numberVectors = vectors[0].dim(0).size;
      type eltType = vectors[0].eltType;
      assert(eltType == real(64));
      createHDF5Dataset(filename, dataset, eltType, (numberVectors, totalCount));

      // Split all states into smaller chunks
      const indexEdges = computeRanges(states, counts);
      const offsets = computeOffsets(indexEdges);
      const numberRanges = indexEdges.domain.dim(1).size - 1;
      // Iterate over chunks
      // TODO: parallelize
      for k in 0 .. numberRanges - 1 {
        // Prepare for kmerge
        var localCounts = [i in indexEdges.dim(0)] indexEdges[i, k + 1] - indexEdges[i, k];
        var maxLocalCount = max reduce localCounts;
        // Copy chunks from different locales to here
        var localStates : [0 .. states.dim(0).size - 1, 0 .. maxLocalCount - 1] uint(64);
        for i in localStates.dim(0) {
          localStates[i, 0 .. localCounts[i] - 1] =
            states[i][indexEdges[i, k] .. indexEdges[i, k + 1] - 1];
        }
        // Run kmerge locally
        const localCount = offsets[k + 1] - offsets[k];
        assert(localCount:int == + reduce localCounts);
        var combinedLocalVectors : [0 .. numberVectors - 1, 0 .. localCount:int - 1] eltType = noinit;
        for ((chunkIndex, indexInChunk), j) in zip(kmergeIndices(localStates, localCounts), 0..) {
          for i in 0 .. numberVectors - 1 {
            combinedLocalVectors[i, j] =
              vectors[chunkIndex][i, indexEdges[chunkIndex, k] + indexInChunk];
          }
        }

        // Save to file
        // TODO: generalize to other eltTypes
        writeHDF5Chunk(filename, dataset, (0:uint, offsets[k]), combinedLocalVectors);
      }
    }
    */

    /*
    module LoadAndDistribute {
      use Time;
      use basis only OnePerLocale;
      use states only hash64_01;
      use wrapper;

      import MatVec.IO.kLoadChunkSize;

      proc splitIntoChunks(lower : int, upper : int, chunkSize : int) {
        assert(upper >= lower);
        assert(chunkSize >= 1);
        const numChunks = (upper - lower + (chunkSize - 1)) / chunkSize;
        var chunks : [0 .. numChunks - 1] (int, int); // = noinit;
        var offset = lower;
        for i in 0 .. numChunks - 1 {
          const nextOffset = min(offset + chunkSize, upper);
          chunks[i] = (offset, nextOffset);
          offset = nextOffset;
        }
        return chunks;
      }

      proc histogramFromChunk(const ref array : [] uint(64), ref timer : Timer) {
        timer.start();
        var histogram : [LocaleSpace] int = 0;
        forall x in array with (+ reduce histogram) {
          const hash = (hash64_01(x) % numLocales:uint):int;
          histogram[hash] += 1;
        }
        timer.stop();
        return histogram;
      }

      private proc checkShapes(totalNumberStates : int, filename : string,
                               datasets : (string, nothing)) {
        const shape = datasetShape(filename, datasets[0]);
        assert(shape.size == 1);
        assert(shape[0] == totalNumberStates);
        return 0;
      }
      private proc checkShapes(totalNumberStates : int, filename : string,
                               datasets : (string, string)) {
        checkShapes(totalNumberStates, filename, (datasets[0], none));
        const shape = datasetShape(filename, datasets[1]);
        assert(shape.size == 2);
        writeln(filename, "/", datasets[1], " has shape ", shape);
        assert(shape[1] == totalNumberStates);
        return shape[0];
      }

      enum WhatToConsider {RepresentativesOnly, ArrayOnly};

      private inline proc determineWhatToLoad(datasets : (string, string)) param : WhatToConsider
      { return WhatToConsider.ArrayOnly; }

      private inline proc determineWhatToLoad(datasets : (string, nothing)) param : WhatToConsider
      { return WhatToConsider.RepresentativesOnly; }

      private proc makeGlobalBuffer(param whatToLoad : WhatToConsider, type eltType,
                                    numberVectors : int, maxCountPerLocale : int) {
        if whatToLoad == WhatToConsider.RepresentativesOnly {
          var globalBuffer : [OnePerLocale] [0 ..# maxCountPerLocale] eltType;
          return globalBuffer;
        }
        else {
          var globalBuffer : [OnePerLocale] [0 ..# numberVectors, 0 ..# maxCountPerLocale] eltType;
          return globalBuffer;
        }
      }

      private proc loadChunks(filename : string, datasets : (string, nothing),
                              type eltType, numberVectors : int, lo : int, hi : int) {
        return (readHDF5Chunk(filename, datasets[0], eltType, (lo,), (hi - lo,)),);
      }
      private proc loadChunks(filename : string, datasets : (string, string),
                              type eltType, numberVectors : int, lo : int, hi : int) {
        const (statesChunk,) =
          loadChunks(filename, (datasets[0], none), uint(64), numberVectors, lo, hi);
        const arrayChunk = readHDF5Chunk(
          filename, datasets[1], eltType, (0, lo), (numberVectors, hi - lo));
        return (statesChunk, arrayChunk);
      }

      private inline proc distributeChunks(ref buffer, ref offset : int,
                                           const ref mask, const ref states) {
        for i in mask.domain {
          if (mask[i]) {
              buffer[here.id][offset] = states[i];
              offset += 1;
          }
        }
      }
      private inline proc distributeChunks(ref buffer, ref offset : int,
                                           const ref mask, const ref states, const ref vectors) {
        const numberVectors = vectors.dim(0).size;
        for i in mask.domain {
          if (mask[i]) {
            for k in 0 ..# numberVectors {
              buffer[here.id][k, offset] = vectors[k, i];
            }
            offset += 1;
          }
        }
      }

      proc mainLoop(filename : string, datasets,
                    type eltType, const ref countPerLocale : [LocaleSpace] int) {
        // Timers
        var _totalTimer : Timer = new Timer();
        _totalTimer.start();
        var _initGlobalBufferTimer = new Timer();
        var _readTimer : [OnePerLocale] Timer = new Timer();
        var _processTimer1 : [OnePerLocale] Timer = new Timer();
        var _copyTimer : [OnePerLocale] Timer = new Timer();

        // Shape info
        const totalNumberStates = + reduce countPerLocale;
        const maxCountPerLocale = max reduce countPerLocale;
        const numberVectors = checkShapes(totalNumberStates, filename, datasets);
        _initGlobalBufferTimer.start();
        var globalBuffer = makeGlobalBuffer(
          determineWhatToLoad(datasets), eltType, numberVectors, maxCountPerLocale);
        _initGlobalBufferTimer.stop();

        coforall loc in Locales /*with (in filename, in datasets)*/
            do on loc {
          const workItems = splitIntoChunks(0, totalNumberStates, kLoadChunkSize);
          var globalOffset : int = 0;
          for (lo, hi) in workItems {
            const size = hi - lo;
            _readTimer[loc.id].start();
            // TODO: why do we need this?
            const localFilename = filename;
            const localDatasets = datasets;
            const chunks = loadChunks(localFilename, localDatasets, eltType, numberVectors, lo, hi);
            _readTimer[loc.id].stop();

            _processTimer1[loc.id].start();
            const ref statesChunk = chunks[0];
            var mask : [statesChunk.domain] bool = noinit;
            forall i in 0 .. size - 1 {
              mask[i] = (hash64_01(statesChunk[i]) % numLocales:uint):int == loc.id;
            }
            _processTimer1[loc.id].stop();

            _copyTimer[loc.id].start();
            distributeChunks(globalBuffer, globalOffset, mask, (...chunks));
            _copyTimer[loc.id].stop();
          }
        }

        _totalTimer.stop();
        writeln("[Chapel] loadAndDistribute took ", _totalTimer.elapsed());
        writeln("[Chapel]     ", _initGlobalBufferTimer.elapsed(),
          " was spent initializing globalBuffer");
        writeln("[Chapel]     ", _readTimer.elapsed(), " was spent reading HDF5 files");
        writeln("[Chapel]     ", _processTimer1.elapsed(), " was spent computing hashes");
        writeln("[Chapel]     ", _copyTimer.elapsed(), " was spent copying to globalBuffer");
        return globalBuffer;
      }
    } // module LoadAndDistribute

    import MatVec.IO.LoadAndDistribute;

    proc calculateNumberStatesPerLocale(filename : string, dataset : string) {
      var _totalTimer = new Timer();
      var _readingTimer : [OnePerLocale] Timer = new Timer();
      var _histogramTimer : [OnePerLocale] Timer = new Timer();
      _totalTimer.start();

      const shape = datasetShape(filename, dataset);
      assert(shape.rank == 1);
      const totalNumberStates = shape[0];
      var histogram : [0 .. numLocales - 1] int;
      const workItems = LoadAndDistribute.splitIntoChunks(0, totalNumberStates, kLoadChunkSize);
      coforall loc in Locales
          with (/*const*/ in workItems, in filename, in dataset,
                + reduce histogram) do on loc {
        for ((lo, hi), chunkIndex) in zip(workItems, 0..) do if (chunkIndex % numLocales == loc.id) {
          _readingTimer[loc.id].start();
          const statesChunk = readHDF5Chunk(filename, dataset, uint(64), (lo,), (hi - lo,));
          _readingTimer[loc.id].stop();
          histogram += LoadAndDistribute.histogramFromChunk(statesChunk, _histogramTimer[loc.id]);
        }
      }
      assert(+ reduce histogram == totalNumberStates);

      _totalTimer.stop();
      writeln("[Chapel] calculateNumberStatesPerLocale took ", _totalTimer.elapsed());
      writeln("[Chapel]     ", _readingTimer.elapsed(), " was spent reading HDF5 files");
      writeln("[Chapel]     ", _histogramTimer.elapsed(), " was spent computing histograms");
      return histogram;
    }

    proc loadRepresentatives(filename : string, basisDataset : string,
                             const ref countPerLocale : [LocaleSpace] int) {
      return LoadAndDistribute.mainLoop(filename, (basisDataset, none), uint(64), countPerLocale);
    }
    proc loadRepresentatives(filename : string, basisDataset : string) {
      return loadRepresentatives(filename, basisDataset,
          calculateNumberStatesPerLocale(filename, basisDataset));
    }

    proc loadVectors(filename : string, basisDataset : string, arrayDataset : string,
                     type eltType, const ref countPerLocale : [LocaleSpace] int) {
      return LoadAndDistribute.mainLoop(filename, (basisDataset, arrayDataset), eltType,
                                        countPerLocale);
    }
    proc loadVectors(filename : string, basisDataset : string, arrayDataset : string,
                     type eltType) {
      return loadVectors(filename, basisDataset, arrayDataset, eltType,
          calculateNumberStatesPerLocale(filename, basisDataset));
    }
    */
  } // module IO
} // module MatVec
