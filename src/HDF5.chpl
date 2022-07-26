module HDF5 {
  use CTypes;
  use BlockDist;

  use LatticeSymmetries.FFI;

  proc datasetRank(filename : string, dataset : string) {
    return ls_hs_hdf5_get_dataset_rank(
      filename.localize().c_str(), dataset.localize().c_str()):int;
  }

  /* Get shape of a dataset in a HDF5 file.
   *
   * :arg filename: Path to HDF5 file.
   * :arg dataset:  Path to dataset within HDF5 file.
   * 
   * :returns: shape of the dataset as a one-dimensional array.
   * :rtype: [] int
   */
  proc datasetShape(filename : string, dataset : string) {
    const rank = datasetRank(filename, dataset);
    var c_shape : [0 ..# rank] uint(64);
    ls_hs_hdf5_get_dataset_shape(filename.localize().c_str(), dataset.localize().c_str(),
      c_ptrTo(c_shape));
    return c_shape:int; 
  }

  private inline proc _makeDomain(shape : 1 * int) : domain(1) { return {0 ..# shape[0]}; }
  private inline proc _makeDomain(shape : 2 * int) : domain(2) { return {0 ..# shape[0], 0 ..# shape[1]}; }

  proc readDataset(filename : string, dataset : string, type eltType, param rank : int)
      where rank == 1 {
    const shape = datasetShape(filename, dataset);
    return readDatasetChunk(filename, dataset, eltType, (0,), (shape[0],));
  }
  proc readDataset(filename : string, dataset : string, type eltType, param rank : int)
      where rank == 2 {
    const shape = datasetShape(filename, dataset);
    return readDatasetChunk(filename, dataset, eltType, (0, 0), (shape[0], shape[1]));
  }

  /* Read part of a dataset from a HDF5 file.
   *
   * :arg filename: Path to HDF5 file.
   * :arg dataset:  Path to dataset within HDF5 file.
   * :arg eltType:  Array datatype.
   * :arg offset:   A tuple of offsets along each dimension.
   * :arg shape:    Array shape.
   *
   * :returns: part of the dataset which is read from file.
   * :rtype: [] eltType
   */
  proc readDatasetChunk(filename : string, dataset : string, type eltType, offset, shape) {
    const dom = _makeDomain(shape);
    var array : [dom] eltType;
    readDatasetChunk(filename, dataset, offset, array);
    return array;
  }

  /* Read part of a dataset from a HDF5 file. This function modifies `array` inplace.
   * 
   */
  proc readDatasetChunk(filename : string, dataset : string, offset, ref array : [] ?eltType)
      where isTuple(offset) && offset.size == array.rank {
    const rank = offset.size;
    var c_offset : [0 .. rank - 1] uint(64) = noinit;
    var c_shape : [0 .. rank - 1] uint(64) = noinit;
    for i in 0 .. rank - 1 {
      c_offset[i] = offset[i]:uint;
      c_shape[i] = array.dim(i).size:uint;
    }
    var args = (filename.localize().c_str(), dataset.localize().c_str(), rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    if (eltType == uint(64)) {
      ls_hs_hdf5_read_chunk_u64((...args));
    }
    else if (eltType == real(64)) {
      ls_hs_hdf5_read_chunk_f64(filename.localize().c_str(), dataset.localize().c_str(),
        rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    }
    else {
      halt("readDatasetChunk does not support " + eltType:string);
    }
  }

  /* Create an HDF5 dataset of given shape and data type.
   *
   */
  proc createHDF5Dataset(filename : string, dataset : string, type eltType, shape) {
    assert(filename.locale == here && dataset.locale == here);
    var c_shape : [0 .. shape.size - 1] uint(64) = noinit;
    for i in 0 .. shape.size - 1 { c_shape[i] = shape[i]:uint; }
    if (eltType == uint(64)) {
      ls_hs_hdf5_create_dataset_u64(filename.localize().c_str(), dataset.localize().c_str(),
        c_shape.size:c_uint, c_ptrTo(c_shape));
    }
    else if (eltType == real(64)) {
      ls_hs_hdf5_create_dataset_f64(filename.localize().c_str(), dataset.localize().c_str(),
        c_shape.size:c_uint, c_ptrTo(c_shape));
    }
    else {
      assert(false);
    }
  }

  /* Write array to a HDF5 dataset.
   */
  proc writeHDF5Chunk(filename : string, dataset : string, offset, array : [?D] ?eltType)
      where isTuple(offset) && offset.size == D.rank {
    assert(D.rank == offset.size);
    var c_offset : [0 .. D.rank - 1] uint(64) = noinit;
    var c_shape  : [0 .. D.rank - 1] uint(64) = noinit;
    for i in c_offset.domain { c_offset[i] = offset[i]:uint; }
    for i in c_shape.domain { c_shape[i] = array.dim(i).size:uint; }

    if (eltType == uint(64)) {
      ls_hs_hdf5_write_chunk_u64(filename.localize().c_str(), dataset.localize().c_str(),
        D.rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    } else if (eltType == real(64)) {
      ls_hs_hdf5_write_chunk_f64(filename.localize().c_str(), dataset.localize().c_str(),
        D.rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    }
    else {
      assert(false);
    }
  }

  proc readBasisStatesAsBlocks(filename : string, dataset : string) {
    const shape = datasetShape(filename, dataset);
    if shape.size != 1 then
      halt("expected '" + dataset + "' to be one-dimensional, but it has shape " + shape:string);
    const totalNumberStates = shape[0];
  
    const box = {0 ..# totalNumberStates};
    const dom = box dmapped Block(box, Locales);
    var states : [dom] uint(64);
    coforall loc in Locales do on loc {
      const indices = states.localSubdomain();
      // logDebug("my subdomain: " + indices:string);
      readDatasetChunk(filename, dataset, (indices.low,), states[indices]);
    }
    return states;
  }

  proc readVectorsAsBlocks(filename : string, dataset : string, type eltType = real(64)) {
    const shape = datasetShape(filename, dataset);
    if shape.size != 2 then
      halt("expected '" + dataset + "' to be two-dimensional, but it has shape " + shape:string);
    const numberVectors = shape[0];
    const totalNumberStates = shape[1];
  
    const boundingBox = {0 ..# numberVectors, 0 ..# totalNumberStates};
    const targetLocales = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
    const dom = boundingBox dmapped Block(boundingBox, targetLocales);
    var vectors : [dom] eltType;
    coforall loc in Locales do on loc {
      const indices = vectors.localSubdomain();
      // logDebug("my subdomain: " + indices:string);
      readDatasetChunk(filename, dataset, indices.low, vectors[indices]);
    }
    return vectors;
  }

}
