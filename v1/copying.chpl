use Time;
use Random;

config const N = 100000;
config const M = 8;
config const iterations = 100;
config param algorithm = 0;

proc distributeBasedOnHash(const ref array : [?D] ?eltType,
                            const ref hashes : [] int,
                            ref timer : Timer) where (D.rank == 2) {
  const numVectors = D.dim(0).size;
  const size = D.dim(1).size;
  var buffer : [0 .. numLocales - 1, 0 .. size - 1, 0 .. numVectors - 1] eltType = noinit;
  var offsets : [LocaleSpace] int;
  for i in 0 .. size - 1 {
    const hash = hashes[i];
    timer.start();
    if (algorithm == 0) {
      for k in 0 .. numVectors - 1 {
        buffer[hash, offsets[hash], k] = array[k, i];
      }
    }
    else if (algorithm == 1) {
      forall k in 0 .. numVectors - 1 {
        buffer[hash, offsets[hash], k] = array[k, i];
      }
    }
    else if (algorithm == 2) {
      const off = offsets[hash];
      buffer[hash, off..off, ..] = array[.., i..i];
    }
    else {
      assert(false);
    }
    timer.stop();
    offsets[hash] += 1;
  }
  return (buffer, offsets);
}

proc main() {
  var array : [0 ..# M, 0 ..# N] real;
  var hashes : [0 ..# N] int;
  fillRandom(hashes);
  hashes %= numLocales;

  var timer = new Timer();
  var (buffer, offsets) = distributeBasedOnHash(array, hashes, timer);
  writeln("Took ", timer.elapsed());
}
