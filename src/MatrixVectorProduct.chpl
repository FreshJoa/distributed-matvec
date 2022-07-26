module MatrixVectorProduct {
  use FFI;
  use ForeignTypes;

  use CTypes;
  use RangeChunk;
  use DynamicIters;

  private proc localCompress(batchSize : int, numberTerms : int, betas : c_ptr(uint(64)),
                             coeffs : c_ptr(?eltType), offsets : c_ptr(int)) {
    var offset : int = 0;
    var i : int = 0;

    offsets[0] = 0;
    for batchIdx in 0 ..# batchSize {
      for termIdx in 0 ..# numberTerms {
        if (coeffs[i] != 0) {
          if (i > offset) {
            coeffs[offset] = coeffs[i];
            betas[offset] = betas[i];
          }
          offset += 1;
        }
        i += 1;
      }
      offsets[batchIdx + 1] = offset;
    }
  }

  // proc localEvaluateWavefunction(const ref alphas : [] uint, const ref state_vector : [] ?eltType,
  //                                ref coeffs : [] ?eltType2, kernels : c_ptr(ls_hs_basis_kernels))
  //     where eltType == eltType2 {
  //   assert(alphas.size == coeffs.size);
  //   const batchSize = alphas.size;
  //   ls_hs_evaluate_wavefunction_via_statevector(
  //     kernels,
  //     batchSize,
  //     c_const_ptrTo(alphas), 1,
  //     c_const_ptrTo(state_vector):c_void_ptr, c_sizeof(eltType),
  //     c_ptrTo(coeffs):c_void_ptr
  //   );
  // }

  record LocalMatVecWorkspace {
    var batchSize : int;
    var numberTerms : int;
    var offDiagDomain : domain(1);
    var spins : [offDiagDomain] uint(64);
    var coeffs : [offDiagDomain] complex(128);
    var offsets : [0 ..# batchSize + 1] int;
    var other_spins : [offDiagDomain] uint(64);
    var other_coeffs : [offDiagDomain] complex(128);
    var norms : [offDiagDomain] real(64);
    var indices : [offDiagDomain] int;

    proc init(batchSize, numberTerms) {
      this.batchSize = batchSize;
      this.numberTerms = numberTerms;
      this.offDiagDomain = {0 ..# (batchSize * numberTerms)};
    }
  }

  private proc computeOffDiag(const ref matrix : Operator,
                              batchSize : int,
                              numberTerms : int,
                              in_spins : c_ptr(uint(64)),
                              out_spins : c_ptr(uint(64)),
                              out_coeffs : c_ptr(complex(128)),
                              out_offsets : c_ptr(int),
                              temp_spins : c_ptr(uint(64)),
                              temp_coeffs : c_ptr(complex(128)),
                              temp_norms : c_ptr(real(64))) {
    if matrix.basis.requiresProjection() {
      ls_hs_operator_apply_off_diag_kernel(matrix.payload, batchSize, in_spins, 1,
        temp_spins, 1, temp_coeffs);
      localCompress(batchSize, numberTerms, temp_spins, temp_coeffs, out_offsets);
      const totalSize = out_offsets[batchSize];
      ls_hs_state_info(matrix.basis.payload, totalSize, temp_spins, 1,
        out_spins, 1, out_coeffs, temp_norms);
      foreach i in 0 ..# totalSize {
        out_coeffs[i] *= temp_coeffs[i] * temp_norms[i];
      }
      ls_hs_state_info(matrix.basis.payload, batchSize, in_spins, 1,
        temp_spins, 1, temp_coeffs, temp_norms);
      for i in 0 ..# batchSize {
        foreach k in out_offsets[i] .. out_offsets[i + 1] - 1 {
          out_coeffs[k] /= temp_norms[i];
        }
      }
    }
    else {
      ls_hs_operator_apply_off_diag_kernel(matrix.payload, batchSize, in_spins, 1,
        out_spins, 1, out_coeffs);
      localCompress(batchSize, numberTerms, out_spins, out_coeffs, out_offsets);
    }
  }

  private proc localEvaluateWaveFunction(const ref matrix : Operator,
                                         const ref xs : [] ?eltType,
                                         batchSize : int,
                                         in_spins : c_ptr(uint(64)),
                                         out_coeffs : c_ptr(?eltType2),
                                         temp_indices : c_ptr(int)) {
    if !matrix.basis.isStateIndexIdentity() {
      ls_hs_state_index(matrix.basis.payload, batchSize,
        in_spins, 1, temp_indices, 1);
      foreach i in 0 ..# batchSize {
        if (temp_indices[i] >= 0) {
          out_coeffs[i] = xs[temp_indices[i]];
        }
        else {
          out_coeffs[i] = 0;
        }
      }
    }
    else {
      foreach i in 0 ..# batchSize {
        out_coeffs[i] = xs[in_spins[i]:int];
      }
    }
  }

  private proc localMatVecPart(const ref matrix : Operator,
                               ref workspace : LocalMatVecWorkspace,
                               startIndex : int, batchSize : int,
                               representatives,
                               xs,
                               ys : [] ?eltType) {
    // Diagonal part
    ls_hs_operator_apply_diag_kernel(matrix.payload, batchSize,
      c_const_ptrTo(representatives[startIndex]), 1, c_ptrTo(workspace.coeffs));
    for i in 0 ..# batchSize {
      const x = xs[startIndex + i];
      const c = workspace.coeffs[i];
      ys[startIndex + i] = (conjg(c) * x):eltType;
    }

    // Off-diagonal part
    const numberTerms = workspace.numberTerms;
    computeOffDiag(matrix, batchSize, numberTerms,
      c_const_ptrTo(representatives[startIndex]),
      c_ptrTo(workspace.spins),
      c_ptrTo(workspace.coeffs),
      c_ptrTo(workspace.offsets),
      c_ptrTo(workspace.other_spins),
      c_ptrTo(workspace.other_coeffs),
      c_ptrTo(workspace.norms));
    localEvaluateWaveFunction(matrix, xs,
      workspace.offsets[batchSize],
      c_ptrTo(workspace.spins),
      c_ptrTo(workspace.other_coeffs),
      c_ptrTo(workspace.indices));
    for i in 0 ..# batchSize {
      var acc : complex(128) = 0;
      for k in workspace.offsets[i] .. workspace.offsets[i + 1] - 1 {
        const x = workspace.other_coeffs[k];
        const c = workspace.coeffs[k];
        acc += conjg(c) * x;
      }
      ys[startIndex + i] += acc:eltType;
    }
  }

  config const localMatrixVectorNumChunks = 100 * here.maxTaskPar;

  proc localMatVec(const ref matrix : Operator,
                   const ref representatives : [] uint(64),
                   const ref xs : [] ?eltType,
                   ref ys : [] eltType) {
   
    const numChunks = min(representatives.size, localMatrixVectorNumChunks);
    const chunkSize = (representatives.size + numChunks - 1) / numChunks;
      // (representatives.size + localMatVecChunkSize - 1) / localMatVecChunkSize;
    var workspace = new LocalMatVecWorkspace(chunkSize, max(1, matrix.numberOffDiagTerms()));
    const todo : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
    // TODO: are we going to run into trouble that there are too many tasks?
    forall _i in dynamic(0 ..# numChunks) with (in workspace) {
      const ref indices = todo[_i];
      assert(indices.size <= chunkSize);
      localMatVecPart(matrix, workspace, indices.low, indices.size,
        representatives, xs, ys);
    }
  }
  proc localMatVec(const ref matrix : Operator,
                   const ref xs : [] ?eltType,
                   ref ys : [] eltType) {
    const ref basis = matrix.basis;
    const ref representatives = basis.representatives();
    localMatVec(matrix, representatives, xs, ys);
  }
}
