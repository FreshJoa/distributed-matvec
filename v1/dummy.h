// This is a workaround before Chapel starts supporting external libraries in
// multilocale lib builds. `ls_flat_spin_basis` is an opaque struct anyway so we
// don't really loose anything.
typedef struct ls_flat_spin_basis ls_flat_spin_basis;
typedef void **ls_void_ptr_ptr;
