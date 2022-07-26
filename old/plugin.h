#pragma once

#include <complex.h>
#include <stdint.h>

void plugin_init(unsigned number_spins);
void plugin_deinit();

typedef struct plugin_build_basis_options_type {
  unsigned number_spins;
  int hamming_weight;
  int spin_inversion;
  int (*is_representative)(void const *internal_state, uint64_t /*count*/,
                           uint64_t const /*spins*/[], uint8_t /*output*/[]);
  void *internal_state;
} plugin_build_basis_options_type;

unsigned plugin_get_number_spins(plugin_build_basis_options_type const *cxt);
int plugin_get_hamming_weight(plugin_build_basis_options_type const *cxt);
int plugin_get_spin_inversion(plugin_build_basis_options_type const *cxt);
int plugin_is_representative(plugin_build_basis_options_type const *cxt,
                             uint64_t count, uint64_t const spins[],
                             uint8_t output[]);

plugin_build_basis_options_type *
plugin_new_options_example_01(unsigned number_spins);
void plugin_delete_options_example_01(plugin_build_basis_options_type *cxt);

uint64_t plugin_get_max_nonzero_per_row();
uint64_t plugin_apply_operator(uint64_t bits, uint64_t *other_spins,
                               _Complex double *other_coeffs);

uint64_t plugin_get_dimension();
uint64_t const *plugin_get_basis_states();
uint64_t plugin_get_index(uint64_t bits);
void plugin_matvec(double const *x, double *y);
