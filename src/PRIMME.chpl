module PRIMME {
  use AllLocalesBarriers;
  use CPtr;
  use SysCTypes;
  require "primme.h";

  extern const PRIMME_VERSION_MAJOR : int;
  extern const PRIMME_VERSION_MINOR : int;

  extern proc primme_initialize(ref primme : primme_params) : void;
  extern proc primme_free(ref primme : primme_params) : void;
  extern proc primme_display_params(in primme : primme_params) : void;

  extern proc sprimme(evals : c_ptr(real(32)), evecs : c_ptr(real(32)),     resNorms : c_ptr(real(32)), const ref primme : primme_params) : c_int;
  extern proc cprimme(evals : c_ptr(real(32)), evecs : c_ptr(complex(64)), resNorms : c_ptr(real(32)), const ref primme : primme_params) : c_int;
  extern proc dprimme(evals : c_ptr(real(64)), evecs : c_ptr(real(64)), resNorms : c_ptr(real(64)), const ref primme : primme_params) : c_int;
  extern proc zprimme(evals : c_ptr(real(64)), evecs : c_ptr(complex(128)), resNorms : c_ptr(real(64)), const ref primme : primme_params) : c_int;

  extern record JD_projectors {
    var LeftQ : c_int;
    var LeftX : c_int;
    var RightQ : c_int;
    var RightX : c_int;
    var SkewQ : c_int;
    var SkewX : c_int;
  }

  extern record correction_params {
    var precondition : c_int;
    var robustShifts : c_int;
    var maxInnerIterations : c_int;
    var projectors : JD_projectors;
    var convTest : primme_convergencetest;
    var relTolBase : real(64);
  }

  // primme_convergencetest enum
  extern type primme_convergencetest = c_int;
  extern const primme_full_LTolerance :primme_convergencetest;
  extern const primme_decreasing_LTolerance :primme_convergencetest;
  extern const primme_adaptive_ETolerance :primme_convergencetest;
  extern const primme_adaptive :primme_convergencetest;

  // primme_event enum
  extern type primme_event = c_int;
  extern const primme_event_outer_iteration :primme_event;
  extern const primme_event_inner_iteration :primme_event;
  extern const primme_event_restart :primme_event;
  extern const primme_event_reset :primme_event;
  extern const primme_event_converged :primme_event;
  extern const primme_event_locked :primme_event;
  extern const primme_event_message :primme_event;
  extern const primme_event_profile :primme_event;

  // primme_init enum
  extern type primme_init = c_int;
  extern const primme_init_default :primme_init;
  extern const primme_init_krylov :primme_init;
  extern const primme_init_random :primme_init;
  extern const primme_init_user :primme_init;

  // primme_op_datatype enum
  extern type primme_op_datatype = c_int;
  extern const primme_op_default :primme_op_datatype;
  extern const primme_op_half :primme_op_datatype;
  extern const primme_op_float :primme_op_datatype;
  extern const primme_op_double :primme_op_datatype;
  extern const primme_op_quad :primme_op_datatype;
  extern const primme_op_int :primme_op_datatype;

  // primme_orth enum
  extern type primme_orth = c_int;
  extern const primme_orth_default :primme_orth;
  extern const primme_orth_implicit_I :primme_orth;
  extern const primme_orth_explicit_I :primme_orth;

  extern record primme_params {
    var n : int(64);
    var matrixMatvec : c_fn_ptr;
    var matrixMatvec_type : primme_op_datatype;
    var applyPreconditioner : c_fn_ptr;
    var applyPreconditioner_type : primme_op_datatype;
    var massMatrixMatvec : c_fn_ptr;
    var massMatrixMatvec_type : primme_op_datatype;
    var numProcs : c_int;
    var procID : c_int;
    var nLocal : int(64);
    var commInfo : c_void_ptr;
    var globalSumReal : c_fn_ptr;
    var globalSumReal_type : primme_op_datatype;
    var broadcastReal : c_fn_ptr;
    var broadcastReal_type : primme_op_datatype;
    var numEvals : c_int;
    var target : primme_target;
    var numTargetShifts : c_int;
    var targetShifts : c_ptr(real(64));
    var dynamicMethodSwitch : c_int;
    var locking : c_int;
    var initSize : c_int;
    var numOrthoConst : c_int;
    var maxBasisSize : c_int;
    var minRestartSize : c_int;
    var maxBlockSize : c_int;
    var maxMatvecs : int(64);
    var maxOuterIterations : int(64);
    var iseed : c_ptr(int(64));
    var aNorm : real(64);
    var BNorm : real(64);
    var invBNorm : real(64);
    var eps : real(64);
    var orth : primme_orth;
    var internalPrecision : primme_op_datatype;
    var printLevel : c_int;
    var outputFile : c_ptr(_file);
    var matrix : c_void_ptr;
    var preconditioner : c_void_ptr;
    var massMatrix : c_void_ptr;
    var ShiftsForPreconditioner : c_ptr(real(64));
    var initBasisMode : primme_init;
    var ldevecs : int(64);
    var ldOPs : int(64);
    var projectionParams : projection_params;
    var restartingParams : restarting_params;
    var correctionParams : correction_params;
    var stats : primme_stats;
    var convTestFun : c_fn_ptr;
    var convTestFun_type : primme_op_datatype;
    var convtest : c_void_ptr;
    var monitorFun : c_fn_ptr;
    var monitorFun_type : primme_op_datatype;
    var monitor : c_void_ptr;
    var queue : c_void_ptr;
    var profile : c_string;
  }

  // primme_preset_method enum
  extern type primme_preset_method = c_int;
  extern const PRIMME_DEFAULT_METHOD :primme_preset_method;
  extern const PRIMME_DYNAMIC :primme_preset_method;
  extern const PRIMME_DEFAULT_MIN_TIME :primme_preset_method;
  extern const PRIMME_DEFAULT_MIN_MATVECS :primme_preset_method;
  extern const PRIMME_Arnoldi :primme_preset_method;
  extern const PRIMME_GD :primme_preset_method;
  extern const PRIMME_GD_plusK :primme_preset_method;
  extern const PRIMME_GD_Olsen_plusK :primme_preset_method;
  extern const PRIMME_JD_Olsen_plusK :primme_preset_method;
  extern const PRIMME_RQI :primme_preset_method;
  extern const PRIMME_JDQR :primme_preset_method;
  extern const PRIMME_JDQMR :primme_preset_method;
  extern const PRIMME_JDQMR_ETol :primme_preset_method;
  extern const PRIMME_STEEPEST_DESCENT :primme_preset_method;
  extern const PRIMME_LOBPCG_OrthoBasis :primme_preset_method;
  extern const PRIMME_LOBPCG_OrthoBasis_Window :primme_preset_method;

  // primme_projection enum
  extern type primme_projection = c_int;
  extern const primme_proj_default :primme_projection;
  extern const primme_proj_RR :primme_projection;
  extern const primme_proj_harmonic :primme_projection;
  extern const primme_proj_refined :primme_projection;

  extern record primme_stats {
    var numOuterIterations : int(64);
    var numRestarts : int(64);
    var numMatvecs : int(64);
    var numPreconds : int(64);
    var numGlobalSum : int(64);
    var numBroadcast : int(64);
    var volumeGlobalSum : int(64);
    var volumeBroadcast : int(64);
    var flopsDense : real(64);
    var numOrthoInnerProds : real(64);
    var elapsedTime : real(64);
    var timeMatvec : real(64);
    var timePrecond : real(64);
    var timeOrtho : real(64);
    var timeGlobalSum : real(64);
    var timeBroadcast : real(64);
    var timeDense : real(64);
    var estimateMinEVal : real(64);
    var estimateMaxEVal : real(64);
    var estimateLargestSVal : real(64);
    var estimateBNorm : real(64);
    var estimateInvBNorm : real(64);
    var maxConvTol : real(64);
    var estimateResidualError : real(64);
    var lockingIssue : int(64);
  }

  // primme_target enum
  extern type primme_target = c_int;
  extern const primme_smallest :primme_target;
  extern const primme_largest :primme_target;
  extern const primme_closest_geq :primme_target;
  extern const primme_closest_leq :primme_target;
  extern const primme_closest_abs :primme_target;
  extern const primme_largest_abs :primme_target;

  // primme_type enum
  extern type primme_type = c_int;
  extern const primme_int :primme_type;
  extern const primme_double :primme_type;
  extern const primme_pointer :primme_type;
  extern const primme_string :primme_type;

  extern record projection_params {
    var projection : primme_projection;
  }

  extern record restarting_params {
    var maxPrevRetain : c_int;
  }

  // proc primmeMatrixMatvec(x : c_void_ptr, ldx : c_ptr(int(64)), y : c_void_ptr, ldy : c_ptr(int(64)), blockSize : c_int,
  //                         primme : c_ptr(primme_params), ierr : c_ptr(c_int)) {
  // }

  // proc eigs(matrix) {
  //   assert(matrix.isLocal);
  //   
  //   var params : primme_params;
  //   params.n = matrix.dimension;
  //   params.numEvals = 1;
  //   params.eps = 1e-9;
  //   params.target = primme_smallest;

  // }

  inline proc _makeInds(shape: int ...?n) {
    var inds : n * range;
    for i in 0 ..# n {
      inds[i] = 0 ..# shape[i];
    }
    return inds;
  }

  pragma "no copy return"
  proc makeArrayFromPtr(ptr : c_ptr, shape)
      where isTuple(shape) && isHomogeneousTuple(shape) && shape[0].type == int {
    var dom = defaultDist.dsiNewRectangularDom(rank=shape.size,
                                               idxType=shape[0].type,
                                               stridable=false,
                                               inds=_makeInds((...shape)));
    dom._free_when_no_arrs = true;
    var arr = new unmanaged DefaultRectangularArr(eltType=ptr.eltType,
                                                  rank=dom.rank,
                                                  idxType=dom.idxType,
                                                  stridable=dom.stridable,
                                                  dom=dom,
                                                  data=ptr:_ddata(ptr.eltType),
                                                  externFreeFunc=nil,
                                                  externArr=true,
                                                  _borrowed=true);
    dom.add_arr(arr, locking = false);
    return _newArray(arr);
  }

  proc main() {
    writeln("Hello world!");

    var a = [1, 2, 3, 4, 5, 6];
    var b : chpl_external_array =
      chpl_make_external_array_ptr(c_ptrTo(a[0]):c_void_ptr, a.size:uint);
    ref c = makeArrayFromPtr(b.elts : c_ptr(int), (2, 3))[.., 0 ..# 2]; // reshape(makeArrayFromExternArray(b, int), {0 ..# 2, 0 ..# 3});
    c[0, 1] = -1;
    writeln(a);
    writeln(c);

    var params : primme_params;
    primme_initialize(params);
    primme_display_params(params);
    primme_free(params);

    /* Set problem matrix */
    // params.matrixMatvec = LaplacianMatrixMatvec;
                            /* Function that implements the matrix-vector product
                               A*x for solving the problem A*x = l*x */
  
    /* Set problem parameters */
    params.n = 100; /* set problem dimension */
    params.numEvals = 10;   /* Number of wanted eigenpairs */
    params.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
    params.target = primme_smallest;
                            /* Wanted the smallest eigenvalues */

    /* Set preconditioner (optional) */
    // primme.applyPreconditioner = LaplacianApplyPreconditioner;
    // params.correctionParams.precondition = 1;
  }
}
