config const n = 3;

const D = {0..<n, 0..<n} dmapped Block(...);  // distribute across all of our locales/compute nodes
var A: [D] real;

var SD: sparse subdomain(D) = {(1, 3), (2, 5), ...};  // COO   // option CSR
var AS: [SD] real;

var x, y: [0..<n] real;

// y = A*x;

forall a in A ...

forall a in A[i..j, ..] ...

forall r in A[r, ..]

forall r in 0..<n do
  y[r] = ...

//
const RowD = {0..<n} dmapped Block(...);
const MatD = {0..<n, 0..<n} dmapped Block(..., targetLocales=.../*make sure we do a row-wise blocking */);

forall i in RowD {
  writeln("locale ", here.id, " owns ", i);
}

forall i in 0..<n do...

begin foo();

on Locales[loc] {
  on A[i,j] {
    ...

    local {
      ...
    }
    
coforall loc in Locales do
  on loc do
    coforall tid in 0..<nTasks do...

forall r in RowD {
  if (debug) then
    if (MatD[r,0].locale != here) then
      halt("Oops, we didn't align our data structures");

  ...MatD.localAccess[r, c]...
}
  
  ...MatD[r, ...]... // will be local to where task is running

// should work today, but wouldn't scale (b/c dim(0) returns a range, which is non-distributed)
forall r in SD.dim(0) do
  y[r] = + reduce [c in SD.dimIter(dim=1, r)] (AS[r,c] * x[c]);

forall r in SD[ do
  y[r] = + reduce [c in SD.dimIter(dim=1, r)] (AS[r,c] * x[c]);

//                   [c in SD[r, ..]]

[i in foo()] expr;

for[all] i in foo() do expr;

// slightly in the future
forall r in SD.dim(0) do
  y[r] = + reduce (AS[r, ..] * x);


var B, C: [1..n] real;

var D = B * C;

forall (b, c) in zip(B, C) do
  b * c;

proc foo(x: real, y: real) {
  ...
}

foo(B, 1.0);

forall b in B do
  foo(b);

// FUTURE SKETCH:
// y = + reduce(dim=columns) AS * x;
