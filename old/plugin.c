#include "plugin.h"
#include <lattice_symmetries/lattice_symmetries.h>

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct plugin_state_type {
  ls_spin_basis *basis;
  ls_operator *matrix;
} plugin_state_type;

plugin_state_type plugin_state;

static void print_error_message_and_exit(ls_error_code const status) {
  fprintf(stderr, "Error (%i): %s\n", status, ls_error_to_string(status));
  exit(1);
}

plugin_build_basis_options_type *
plugin_new_options_example_01(unsigned const number_spins) {
  ls_enable_logging();
  LATTICE_SYMMETRIES_LOG_DEBUG("Creating plugin_build_basis_options_type for a "
                               "linear chain of %u spins...\n",
                               number_spins);
  int const hamming_weight = number_spins / 2;
  ls_group *group;
  ls_error_code status = ls_create_group(&group, 0, NULL);
  if (status != LS_SUCCESS) {
    goto fail1;
  }

  ls_spin_basis *basis;
  status = ls_create_spin_basis(&basis, group, number_spins, hamming_weight,
                                /*spin_inversion=*/0);
  if (status != LS_SUCCESS) {
    goto fail2;
  }

  LATTICE_SYMMETRIES_LOG_DEBUG(
      "%s\n",
      "Building a list of representatives on C side. This may take a while...");
  status = ls_build(basis);
  if (status != LS_SUCCESS) {
    goto fail3;
  }

  plugin_build_basis_options_type *cxt =
      malloc(sizeof(plugin_build_basis_options_type));
  if (cxt == NULL) {
    status = LS_OUT_OF_MEMORY;
    goto fail4;
  }

  cxt->number_spins = ls_get_number_spins(basis);
  cxt->hamming_weight = ls_get_hamming_weight(basis);
  cxt->spin_inversion = ls_get_spin_inversion(basis);
  cxt->is_representative = ls_is_representative;
  cxt->internal_state = basis;

fail4:
fail3:
  if (status != LS_SUCCESS) {
    ls_destroy_spin_basis(basis);
  }
fail2:
  ls_destroy_group(group);
fail1:
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }

  return cxt;
}

void plugin_delete_options_example_01(plugin_build_basis_options_type *cxt) {
  cxt->number_spins = 0;
  cxt->hamming_weight = 9;
  cxt->spin_inversion = 0;
  cxt->is_representative = NULL;
  ls_destroy_spin_basis(cxt->internal_state);
}

#if 0
void plugin_init(unsigned const number_spins) {
  int const hamming_weight = number_spins / 2;

  ls_group *group;
  ls_error_code status = ls_create_group(&group, 0, NULL);
  if (status != LS_SUCCESS) {
    goto fail1;
  }

  ls_spin_basis *basis;
  status = ls_create_spin_basis(&basis, group, number_spins, hamming_weight,
                                /*spin_inversion=*/0);
  if (status != LS_SUCCESS) {
    goto fail2;
  }
  status = ls_build(basis);
  if (status != LS_SUCCESS) {
    goto fail3;
  }

  // Get number of representatives
  uint64_t number_states;
  status = ls_get_number_states(basis, &number_states);
  if (status != LS_SUCCESS) {
    goto fail4;
  }
  printf("Hilbert space dimension is %zu\n", (size_t)number_states);

  // Heisenberg Hamiltonian
  // clang-format off
  _Complex double const matrix[] = {1.0,  0.0,  0.0, 0.0,
                                    0.0, -1.0,  2.0, 0.0,
                                    0.0,  2.0, -1.0, 0.0,
                                    0.0,  0.0,  0.0, 1.0};
  // clang-format on
  uint16_t(*edges)[2] = malloc(number_spins * sizeof(uint16_t[2]));
  if (edges == NULL) {
    status = LS_OUT_OF_MEMORY;
    goto fail5;
  }
  for (unsigned i = 0U; i < number_spins; ++i) {
    edges[i][0] = i;
    edges[i][1] = (i + 1) % number_spins;
  }
  ls_interaction *term;
  status = ls_create_interaction2(&term, matrix, number_spins, edges);
  if (status != LS_SUCCESS) {
    goto fail6;
  }
  ls_operator *hamiltonian;
  status = ls_create_operator(&hamiltonian, basis, 1,
                              (ls_interaction const *const *)&term);
  if (status != LS_SUCCESS) {
    goto fail7;
  }

fail7:
  ls_destroy_interaction(term);
fail6:
  free(edges);
fail5:
fail4:
fail3:
  if (status != LS_SUCCESS) {
    ls_destroy_spin_basis(basis);
  }
fail2:
  ls_destroy_group(group);
fail1:
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }

  plugin_state = (plugin_state_type){basis, hamiltonian};
}

void plugin_deinit() {
  ls_destroy_spin_basis(plugin_state.basis);
  ls_destroy_operator(plugin_state.matrix);
}
#endif

unsigned plugin_get_number_spins(plugin_build_basis_options_type const *cxt) {
  return cxt->number_spins;
}

int plugin_get_hamming_weight(plugin_build_basis_options_type const *cxt) {
  return cxt->hamming_weight;
}

int plugin_get_spin_inversion(plugin_build_basis_options_type const *cxt) {
  return cxt->spin_inversion;
}

int plugin_is_representative(plugin_build_basis_options_type const *cxt,
                             uint64_t const count, uint64_t const spins[],
                             uint8_t output[]) {
  return cxt->is_representative(cxt->internal_state, count, spins, output);
}

#if 0
uint64_t plugin_get_dimension() {
  uint64_t dim;
  ls_error_code status = ls_get_number_states(plugin_state.basis, &dim);
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }
  return dim;
}

uint64_t plugin_get_max_nonzero_per_row() {
  return ls_operator_max_buffer_size(plugin_state.matrix);
}

uint64_t const *plugin_get_basis_states() {
  ls_states *states;
  ls_get_states(&states, plugin_state.basis);
  uint64_t const *p = ls_states_get_data(states);
  ls_destroy_states(states);
  return p;
}

typedef struct plugin_apply_callback_state_type {
  uint64_t *spins;
  _Complex double *coeffs;
  uint64_t size;
} plugin_apply_callback_state_type;

static ls_error_code plugin_apply_callback(ls_bits512 const *bits,
                                           void const *coeff, void *cxt) {
  plugin_apply_callback_state_type *state = cxt;
  state->spins[state->size] = bits->words[0];
  state->coeffs[state->size] = *((_Complex double const *)coeff);
  ++(state->size);
  return LS_SUCCESS;
}

uint64_t plugin_apply_operator(uint64_t const bits, uint64_t *other_spins,
                               _Complex double *other_coeffs) {
  ls_bits512 spin;
  memset(&spin, 0, sizeof(ls_bits512));
  spin.words[0] = bits;

  plugin_apply_callback_state_type state = {other_spins, other_coeffs, 0};
  ls_error_code status = ls_operator_apply(plugin_state.matrix, &spin,
                                           plugin_apply_callback, &state);
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }
  return state.size;
}

uint64_t plugin_get_index(uint64_t const bits) {
  uint64_t index;
  ls_error_code status = ls_get_index(plugin_state.basis, bits, &index);
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }
  return index;
}

void plugin_matvec(double const *x, double *y) {
  uint64_t size;
  ls_get_number_states(plugin_state.basis, &size);
  ls_error_code status =
      ls_operator_matmat(plugin_state.matrix, LS_FLOAT64, size, 1, x, 1, y, 1);
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }
}
#endif
