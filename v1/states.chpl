use BitOps;
use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use Time;

use wrapper;
use profiling;

/* Hash function which we use to map spin configurations to locale indices.

   Typical usage:
   ```chapel
   var localeIndex = (hash64_01(x) % numLocales:uint):int;
   ```
*/
inline proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

// proc hash32_01(in x: uint(32)): uint(32) {
//   x += 1;
//   x ^= x >> 17;
//   x *= 0xed5ad4bb:uint(32);
//   x ^= x >> 11;
//   x *= 0xac4c1b51:uint(32);
//   x ^= x >> 15;
//   x *= 0x31848bab:uint(32);
//   x ^= x >> 14;
//   return x;
// }

/* Get next integer with the same hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
inline proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

inline proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

inline proc nextState(v: uint(64), isHammingWeightFixed : bool): uint(64) {
  return if isHammingWeightFixed then nextStateFixedHamming(v) else nextStateGeneral(v);
}

/* Fill buffer using either nextStateFixedHamming or nextStateGeneral. Returns
   the number of elements written.
 */
inline proc manyNextState(in v: uint(64), bound: uint(64), buffer: [] uint(64),
                          isHammingWeightFixed : bool) {
  assert(v <= bound);
  for i in buffer.domain {
    buffer[i] = v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if (v == bound) { return i + 1; }
    v = nextState(v, isHammingWeightFixed);
  }
  return buffer.size;
}

inline iter findStatesInRange(in lower: uint(64), upper: uint(64),
                              isHammingWeightFixed : bool,
                              // nextStateFn: func(uint(64), uint(64)),
                              basis: c_ptr(ls_spin_basis)) {
  /* C code will use SIMD, and it's beneficial to call ls_is_representative with
     count > 1. chunkSize controls the size of chunks which will be fed to
     ls_is_representative.
   */
  const chunkSize = 10240;
  var buffer: [0 .. chunkSize - 1] uint(64);
  var flags: [0 .. chunkSize - 1] uint(8);
  while (true) {
    const written = manyNextState(lower, upper, buffer, isHammingWeightFixed);
    ls_is_representative(
        basis, written:uint, c_ptrTo(buffer), c_ptrTo(flags));
    for i in 0 ..# written {
      if (flags[i]:bool) { yield buffer[i]; }
    }
    if (buffer[written - 1] == upper) { break; }
    lower = nextState(buffer[written - 1], isHammingWeightFixed);
  }
}

inline proc testBit(bits: uint(64), i: int): bool { return ((bits >> i) & 1):bool; }
inline proc clearBit(bits: uint(64), i: int) { return bits & ~(1:uint(64) << i); }
inline proc setBit(bits: uint(64), i: int) { return bits | (1:uint(64) << i); }

/* Get the closest to `x` integer with Hamming weight `hammingWeight`.
 */
inline proc closestWithFixedHamming(in x: uint(64), hammingWeight: uint): uint(64) {
  assert(hammingWeight <= 64);
  var weight = popcount(x);
  if (weight > hammingWeight) {
      // Keep clearing lowest bits until we reach the desired Hamming weight.
      var i = 0;
      while (weight > hammingWeight) {
        assert(i < 64);
        if (testBit(x, i)) {
            x = clearBit(x, i);
            weight -= 1;
        }
        i += 1;
      }
      // NOTE: why??
      // var maxValue = hammingWeight == 0 ? (0:uint(64)) : (~(0:uint(64)) << (64 - hammingWeight));
      // if (x < max_value) { x = next_state<true>(x); }
  }
  else if (weight < hammingWeight) {
    // Keep setting lowest bits until we reach the desired Hamming weight.
    var i = 0;
    while (weight < hammingWeight) {
      assert(i < 64);
      if (!testBit(x, i)) {
          x = setBit(x, i);
          weight += 1;
      }
      i += 1;
    }
  }
  return x;
}

var _splitIntoRangesTime = new MeasurementTable("splitIntoRanges");

/* Splits `[current, bound]` into chunks
   `{[current, upper_1], [lower_2, upper_2], ..., [lower_N, bound]}`
   such that `upper_i - lower_i` is approximately `chunkSize`. Note that care
   is taken to ensure that each `lower_i` and `upper_i` has the right Hamming
   weight when `isHammingWeightFixed` is `true`.
 */
proc splitIntoRanges(in current: uint(64), bound: uint(64), chunkSize: uint(64),
                     isHammingWeightFixed: bool) {
  var __timer = getTimerFor(_splitIntoRangesTime);
  var ranges: list((uint(64), uint(64)));
  const hammingWeight = popcount(current);
  while (true) {
      if (bound - current <= chunkSize) {
          ranges.append((current, bound));
          break;
      }
      var next = current + chunkSize;
      if isHammingWeightFixed { next = closestWithFixedHamming(next, hammingWeight); }

      assert(next >= current);
      if (next >= bound) {
          ranges.append((current, bound));
          break;
      }
      ranges.append((current, next));
      current = if isHammingWeightFixed then nextStateFixedHamming(next)
                                        else nextStateGeneral(next);
  }
  var distRanges = newCyclicArr({0.. ranges.size - 1}, (uint(64), uint(64)));
  distRanges = ranges;
  return distRanges;
}
