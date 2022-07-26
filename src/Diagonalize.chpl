module Diagonalize {
  use PRIMME;
  use ApplyOperator;

  // A buffer located on locale 0 to help with the broadcast                      
  var tmpBuffDom = {0..0:c_int};
  var tmpBuff: [tmpBuffDom] real;

  proc broadcastReal(buffer: c_ptr(real), count: c_ptr(c_int)) {
    const n = count.deref(),
          inds = 0..<n;

    if here.id == 0 {
      // grow the temp buff if it's not big enough                                
      if n > tmpBuffDom.size then
        tmpBuffDom = {inds};

      // copy locale 0's data into the buffer                                     
      forall i in inds do
        tmpBuff[i] = buffer[i];
    }

    // wait until locale 0's got tmpBuff set up before proceeding                 
    allLocalesBarrier.barrier();

    // Locale 0 already has the data so doesn't need to do anything               
    if (here.id != 0) then
      forall i in inds do
        buffer[i] = tmpBuff[i];
  }

  export proc ls_chpl_broadcast_real(sendBuf : c_void_ptr, count : c_ptr(c_int),
                                     primme : c_ptr(primme_params), ierr : c_ptr(c_int)) {
    if primme.deref().broadcastReal_type != primme_op_double then
      halt("broadcastReal is implemented for double precision only");
    broadcastReal(sendBuf : c_ptr(real), count);
    ierr.deref() = 0;
  }

  // A buffer of atomics on locale 0 for computing the reduction                  
  var atomicBuffDom = {0..0:c_int};
  var atomicBuff: [atomicBuffDom] atomic real;

  proc globalSumReal(sendBuf: c_ptr(real), recvBuf: c_ptr(real),
                     count: c_ptr(c_int)) {
    const n = count.deref(),
          inds = 0..<n;

    // grow the temp buff if it's not big enough                                  
    if here.id == 0 then
      if n > atomicBuffDom.size then
        atomicBuffDom = {inds};

    // Make sure locale 0 has had the chance to resize before proceeding          
    allLocalesBarrier.barrier();

    // have all locales atomically add their results to the atomicBuff            
    forall i in inds do
      atomicBuff[i].add(sendBuf[i]);

    // Make sure all locales have accumulated their contributions                 
    allLocalesBarrier.barrier();

    // Have each locale copy the results out into its buffer                      
    forall i in inds do
      recvBuf[i] = atomicBuff[i].read();
  }

  export proc ls_chpl_global_sum_real(sendBuf : c_void_ptr, recvBuf : c_void_ptr, count : c_ptr(c_int),
                                      primme : c_ptr(primme_params), ierr : c_ptr(c_int)) {
    if primme.deref().globalSumReal_type != primme_op_double then
      halt("globalSum is implemented for double precision only");
    globalSumReal(sendBuf : c_ptr(real), recvBuf : c_ptr(real), count);
    ierr.deref() = 0;
  }

  inline proc borrowOperator(primme : c_ptr(primme_params)) {
    const p = primme.deref().matrix : c_ptr(ls_hs_operator);
    return new Operator(p, owning=false);
  }

  export proc ls_chpl_primme_matvec(_x : c_void_ptr, _ldx : c_ptr(int(64)), _y : c_void_ptr, _ldy : c_ptr(int(64)),
                                    _blockSize : c_ptr(c_int), primme : c_ptr(primme_params), _ierr : c_ptr(c_int)) {
    const blockSize = _blockSize.deref():int;
    if blockSize != 1 then
      halt("currently only blockSize=1 is supported");
    const ldx = _ldx.deref();
    const ldy = _ldy.deref();
    const n = primme.deref().n;
    assert(ldx >= n);
    assert(ldy >= n);

    // const precision = primme.deref().matrixMatvec_type;
    // const dtype = 
    type eltType = real;

    ref X = makeArrayFromPtr(_x : c_ptr(eltType), (blockSize, ldx))[.., 0 ..# n];
    ref Y = makeArrayFromPtr(_y : c_ptr(eltType), (blockSize, ldy))[.., 0 ..# n];
    var matrix = borrowOperator(primme);
    localMatVec(matrix, X, Y);

    _ierr.deref() = 0;
  }
}
