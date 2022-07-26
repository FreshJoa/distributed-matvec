use CPtr;
use SysCTypes;
use profiling;

require "lattice_symmetries/lattice_symmetries.h";
extern type ls_spin_basis;
extern type ls_operator;
extern type ls_flat_spin_basis;
extern type ls_error_code;

extern proc ls_enable_logging();
extern proc ls_get_number_spins(basis: c_ptr(ls_spin_basis)): c_uint;
extern proc ls_get_hamming_weight(basis: c_ptr(ls_spin_basis)): c_int;
extern proc ls_get_spin_inversion(basis: c_ptr(ls_spin_basis)): c_int;
extern proc ls_is_representative(basis: c_ptr(ls_spin_basis), count: uint(64),
                                 bits: c_ptr(uint(64)), is_repr: c_ptr(uint(8)));
extern proc ls_operator_max_buffer_size(op : c_ptr(ls_operator)): uint(64);
extern proc ls_binary_search(data : c_ptr(uint(64)), size : uint(64), key : uint(64)): uint(64);
// extern proc ls_convert_to_flat_spin_basis(ptr: c_ptr(c_ptr(ls_flat_spin_basis)),
//                                           basis: c_ptr(ls_spin_basis)): ls_error_code;
// extern proc ls_destroy_flat_spin_basis(ptr: c_ptr(ls_flat_spin_basis));
// extern proc ls_get_buffer_size_for_flat_spin_basis(basis: c_ptr(ls_flat_spin_basis)): uint(64);
// extern proc ls_serialize_flat_spin_basis(basis: c_ptr(ls_flat_spin_basis),
//                                          buffer: [] c_char, size: uint(64)): ls_error_code;
// extern proc ls_deserialize_flat_spin_basis(ptr: c_ptr(c_ptr(ls_flat_spin_basis)),
//                                            buffer: [] c_char, size: uint(64)): ls_error_code ;
// extern proc ls_copy_spin_basis(basis: c_ptr(ls_spin_basis)): c_ptr(ls_spin_basis);

require "lattice_symmetries_haskell.h";
extern record ls_hs_spin_basis_v1 {
  var payload: c_ptr(ls_spin_basis);
  var context: c_void_ptr;
}

extern record ls_hs_operator_v1 {
  var payload: c_ptr(ls_operator);
  var context: c_void_ptr;
}

extern proc ls_hs_init();
extern proc ls_hs_exit();
extern proc ls_hs_basis_and_hamiltonian_from_yaml(path: c_string,
    basis: c_ptr(ls_hs_spin_basis_v1), hamiltonian: c_ptr(ls_hs_operator_v1));
extern proc ls_hs_destroy_spin_basis(basis: c_ptr(ls_hs_spin_basis_v1));
extern proc ls_hs_destroy_operator(basis: c_ptr(ls_hs_operator_v1));
extern proc ls_hs_operator_apply(op : c_ptr(ls_hs_operator_v1),
    count : uint(64), spins : c_ptr(uint(64)), offsets : c_ptr(uint(64)),
    out_spins : c_ptr(uint(64)), out_coeffs : c_ptr(complex(128))) : c_int;
extern proc ls_hs_hdf5_get_dataset_rank(path: c_string, dataset: c_string):c_uint;
extern proc ls_hs_hdf5_get_dataset_shape(path: c_string, dataset: c_string,
                                         shape: c_ptr(uint(64)));
extern proc ls_hs_hdf5_create_dataset_u64(path: c_string,
    dataset: c_string, dim: c_uint, shape: c_ptr(uint(64)));
extern proc ls_hs_hdf5_create_dataset_f64(path: c_string,
    dataset: c_string, dim: c_uint, shape: c_ptr(uint(64)));
extern proc ls_hs_hdf5_write_chunk_u64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(uint(64)));
extern proc ls_hs_hdf5_write_chunk_f64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(real(64)));
extern proc ls_hs_hdf5_read_chunk_u64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(uint(64)));
extern proc ls_hs_hdf5_read_chunk_f64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(real(64)));

proc initRuntime(withLogging : bool) {
  for loc in Locales do on loc {
    if (withLogging) { ls_enable_logging(); }
    // writeln("[Chapel] Calling ls_hs_init on Locales[", loc.id, "] ...");
    ls_hs_init();
    // writeln("[Chapel] ls_hs_init on Locales[", loc.id, "] is done");
  }
}

proc deinitRuntime() {
  for loc in Locales do on loc {
    // writeln("[Chapel] Calling ls_hs_exit on Locales[", loc.id, "] ...");
    ls_hs_exit();
    // writeln("[Chapel] ls_hs_exit on Locales[", loc.id, "] is done");
  }
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
  assert(filename.locale == here && dataset.locale == here);
  const rank = ls_hs_hdf5_get_dataset_rank(
      filename.localize().c_str(), dataset.localize().c_str()):int;
  var c_shape : [0 .. rank - 1] uint(64);
  ls_hs_hdf5_get_dataset_shape(filename.localize().c_str(), dataset.localize().c_str(),
      c_ptrTo(c_shape));
  return [i in c_shape.domain] c_shape[i]:int; 
}

proc _makeDomain(shape : 1 * int) : domain(1) { return {0 .. shape[0]:int - 1}; }
proc _makeDomain(shape : 2 * int) : domain(2) {
  return {0 .. shape[0]:int - 1, 0 .. shape[1]:int - 1};
}

var _readHDF5ChunkTime = new MeasurementTable("readHDF5Chunk");

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
proc readHDF5Chunk(filename : string, dataset : string, type eltType, offset, shape) {
  const dom = _makeDomain(shape);
  var array : [dom] eltType;
  readHDF5Chunk(filename, dataset, offset, array);
  return array;
}

/* Read part of a dataset from a HDF5 file. This function modifies `array` inplace.
 * 
 */
proc readHDF5Chunk(filename : string, dataset : string, offset, array : [] ?eltType)
    where isTuple(offset) && offset.size == array.rank {
  var __timer = getTimerFor(_readHDF5ChunkTime);
  const rank = offset.size;
  var c_offset : [0 .. rank - 1] uint(64) = noinit;
  var c_shape : [0 .. rank - 1] uint(64) = noinit;
  for i in 0 .. rank - 1 {
    c_offset[i] = offset[i]:uint;
    c_shape[i] = array.dim(i).size:uint;
  }
  if (eltType == uint(64)) {
    ls_hs_hdf5_read_chunk_u64(filename.localize().c_str(), dataset.localize().c_str(),
      rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
  }
  else if (eltType == real(64)) {
    ls_hs_hdf5_read_chunk_f64(filename.localize().c_str(), dataset.localize().c_str(),
      rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
  }
  else {
    assert(false);
  }
  // return array;
}

var _createHDF5DatasetTime = new MeasurementTable("createHDF5Dataset");

/* Create an HDF5 dataset of given shape and data type.
 *
 */
proc createHDF5Dataset(filename : string, dataset : string, type eltType, shape) {
  assert(filename.locale == here && dataset.locale == here);
  var __timer = getTimerFor(_createHDF5DatasetTime);
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

var _writeHDF5ChunkTime = new MeasurementTable("writeHDF5Chunk");

/* Write array to a HDF5 dataset.
 */
proc writeHDF5Chunk(filename : string, dataset : string, offset, array : [?D] ?eltType)
    where isTuple(offset) && offset.size == D.rank {
  assert(D.rank == offset.size);
  var __timer = getTimerFor(_writeHDF5ChunkTime);
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
