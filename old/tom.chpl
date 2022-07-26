config const n = 3;


const D = {1..n, 1..n};

/*
var A: [D] real;

var A: [1..n, 1..n] real;

writeln(A);

forall a in A do
  a += 0.1;

writeln(A);

forall (i,j) in A.domain do
  A[i,j] += 0.1;

writeln(A);
*/


proc B(i: int, j: int, type eltType = real) {
  return (i*1000 + j:eltType/10): eltType;
}


writeln(B(4,7));
writeln(B[4,7]);

for (i,j) in D do
  writeln("B", (i,j) , " = ", B(i,j));

