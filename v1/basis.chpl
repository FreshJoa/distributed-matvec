use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use Time;
use VisualDebug;
use RangeChunk;
use Search;
use DynamicIters;

use wrapper;
use states;
use profiling;
use Merge;

config const kNoSymmetries = false;

const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);

class DistributedBasis {
  var _localBases: [OnePerLocale] ls_hs_spin_basis_v1;

  proc init(path: string) {
    this._localBases = [loc in Locales] new ls_hs_spin_basis_v1(nil, nil);
    complete();

    coforall loc in Locales do on loc {
      var dummyHamiltonian = new ls_hs_operator_v1(nil, nil);
      ls_hs_basis_and_hamiltonian_from_yaml(
        path.localize().c_str(), c_ptrTo(this._localBases[loc.id]), c_ptrTo(dummyHamiltonian));
      ls_hs_destroy_operator(c_ptrTo(dummyHamiltonian));
    }
  }

  proc deinit() {
    coforall loc in Locales do on loc {
      ls_hs_destroy_spin_basis(c_ptrTo(this._localBases[loc.id]));
    }
  }

  inline proc rawPtr() { return this._localBases[here.id].payload; }
  inline proc numberSpins() { return ls_get_number_spins(this.rawPtr()); }
  inline proc hammingWeight() { return ls_get_hamming_weight(this.rawPtr()); }
  inline proc spinInversion() { return ls_get_spin_inversion(this.rawPtr()); }
  inline proc isHammingWeightFixed() { return this.hammingWeight() != -1; }

  proc statesBounds() {
    var lower: uint(64);
    var upper: uint(64);
    if (this.isHammingWeightFixed()) {
      lower = ~(0:uint(64)) >> (64 - this.hammingWeight());
      upper = lower << (this.numberSpins() - this.hammingWeight());
    }
    else {
      lower = 0;
      if (this.numberSpins() == 64) { upper = ~(0:uint(64)); }
      else { upper = (1:uint(64) << this.numberSpins()) - 1; }
    }
    if spinInversion() != 0 {
      upper = upper / 2;
    }
    assert(lower <= upper);
    return (lower, upper);
  }
}


inline proc makeShift(numberSpins : int, numberBits : int) : int {
  assert(0 < numberBits && numberBits < 64);
  return if numberBits >= numberSpins then 0 else (numberSpins - numberBits);
}

var _getIndexTime = new MeasurementTable("BasisStates::getIndex");

class BasisStates {
  var _size : int;
  var _representatives : [OnePerLocale] [0 ..# _size] uint(64);
  var _counts : [OnePerLocale] int;
  var _numberSpins : int;
  var _numberBits : int = 16;
  var _shift : int;
  var _ranges : [OnePerLocale] [0 ..# 1 + (1 << _numberBits)] int;

  proc _generateBucketRanges() {
    coforall loc in Locales do on loc {
      ref localRanges = _ranges[loc.id];
      ref localStates = _representatives[loc.id][0 ..# _counts[loc.id]];
      const numberBuckets = 1 << _numberBits;
      var offset = 0;
      for i in 0 ..# numberBuckets {
        localRanges[i] = offset;
        // writeln("ranges[", loc.id, "][", i, "] = ", offset);
        while offset != localStates.size && bucketIndex(localStates[offset]) == i {
          offset += 1;
        }
      }
      localRanges[numberBuckets] = offset;
      assert(offset == localStates.size);
    }
  }

  proc init(numberSpins : int, representatives, counts, numberBits : int) {
    this._size = representatives[0].size;
    this._representatives = representatives;
    this._counts = counts;
    this._numberSpins = numberSpins;
    this._numberBits = min(numberSpins, numberBits);
    this._shift = makeShift(numberSpins, numberBits);
    this.complete();
    _generateBucketRanges();
  }

  inline proc totalNumberStates() : int {
    return + reduce _counts;
  }

  inline proc bucketIndex(x : uint(64)) : int {
    return (x >> _shift):int;
  }

  inline proc getCounts(loc : locale = here) : int {
    return _counts[loc.id];
    // ref r = _representatives[loc.id];
    // return r[0 ..# _counts[loc.id]];
  }

  inline proc getRepresentatives(loc : locale = here) ref {
    return _representatives[loc.id];
    // ref r = _representatives[loc.id][0 ..# _counts[loc.id]];
    // return r;
  }

  proc getIndexFn() {
    record Getter {
      const shift : int;
      const rangesPtr: c_ptr(int);
      const dataPtr: c_ptr(uint(64));

      inline proc this(x : uint(64)) : int {
        if (kNoSymmetries) { assert(dataPtr[x:int] == x); return x:int; }
        const i = (x >> shift):int;
        const b = rangesPtr[i];
        const e = rangesPtr[i + 1];
        const size = (e - b):uint(64);
        const j = ls_binary_search(dataPtr + b, size, x);
        assert(j < size);
        return b + j:int;
      }
    }
    return new Getter(_shift, c_ptrTo(_ranges[here.id]), c_ptrTo(_representatives[here.id]));
  }

  proc getIndex(x : uint(64)) : int {
    // var __time = getTimerFor(_getIndexTime);
    const i = bucketIndex(x);
    const b = _ranges[here.id][i];
    const e = _ranges[here.id][i + 1];
    // The following causes a memory allocation... :(
    // ref r = getRepresentatives(loc)[_ranges[loc.id][i] .. _ranges[loc.id][i + 1] - 1];
    const dataPtr = c_ptrTo(_representatives[here.id][b]);
    const size = (e - b):uint(64);
    var j : uint(64) = ls_binary_search(dataPtr, size, x);
    assert(j < size);
    return b + j:int;
  }
}

// proc numlocs() param where (CHPL_COMM==none) {
//   return 1;
// }
// proc numlocs() {
//   return numLocales;
// }

config const kMakeStatesMinChunkSize = 1024 * 1024;
config const kMakeStatesMaxNumChunks = numLocales * here.maxTaskPar * 5;

var _makeStatesTime = new MeasurementTable("makeStates");
var _makeStatesLoopTime = new MeasurementTable("makeStates::main loop");
var _makeStatesAllocatingTime = new MeasurementTable("makeStates::allocating chunks");
var _makeStatesListToArrayTime = new MeasurementTable("makeStates::list to array");
var _makeStatesRemoteCopiesTime = new MeasurementTable("makeStates::remote copy");

proc makeStates(basis: DistributedBasis, bits : int = 16) {
  var __timer = getTimerFor(_makeStatesTime);
  const (lower, upper) = basis.statesBounds();
  const chunkSize =
    max(kMakeStatesMinChunkSize:uint, (upper - lower) / kMakeStatesMaxNumChunks:uint);
  const ranges = splitIntoRanges(lower, upper, chunkSize, basis.isHammingWeightFixed());
  writeln("[Chapel] makeStates is using chunkSize = ", chunkSize,
          "; constructed ", ranges.size, " ranges");
  var __timer1 = getTimerFor(_makeStatesAllocatingTime);
  const D = {0 .. ranges.size - 1} dmapped Cyclic(startIdx=0);
  var chunks: [D] [LocaleSpace] list(uint(64));
  __timer1.stop();

  // This is the part that should take the most amount of time in this function.
  // Everything else is pure overhead...
  var __timer2 = getTimerFor(_makeStatesLoopTime);
  coforall loc in Locales do on loc {
    var localIndices = chunks.localSubdomain(loc);
    forall i in dynamic(localIndices) {
      // coforall tid in 0 ..# here.maxTaskPar {
      // var taskIndices = localIndices by here.maxTaskPar align (tid * numLocales + loc.id);
      // for i in taskIndices {
      const (lower, upper) = ranges[i];
      ref chunk = chunks[i];
      for x in findStatesInRange(lower, upper, basis.isHammingWeightFixed(),
                                 basis.rawPtr()) {
        const hash = (hash64_01(x) % numLocales:uint):int;
        chunk[hash].append(x);
      }
      // }
    }
  }
  __timer2.stop();

  var __timer3 = getTimerFor(_makeStatesListToArrayTime);
  var maxPerChunk = max reduce ([i in D] (
                      max reduce ([j in LocaleSpace]
                        chunks[i][j].size)));
  var arrChunks: [D] [LocaleSpace] [0 ..# maxPerChunk] uint(64);
  var arrSizes: [D] [LocaleSpace] int;
  forall i in D {
    forall j in LocaleSpace {
      const n = chunks[i][j].size;
      arrChunks[i][j][0 ..# n] = chunks[i][j];
      arrSizes[i][j] = n;
    }
  }
  __timer3.stop();

  var __timer4 = getTimerFor(_makeStatesRemoteCopiesTime);
  var counts: [OnePerLocale] int =
    [j in LocaleSpace] (+ reduce [i in 0..<ranges.size] chunks[i][j].size);
  const maxCount = max reduce counts;
  var states: [OnePerLocale] [0 ..# maxCount] uint(64);
  coforall loc in Locales do on loc {
    var offset: int = 0;
    for i in 0 ..# ranges.size {
      var c = arrSizes[i][loc.id]; // chunks[i][loc.id].size;
      states[loc.id][offset ..# c] = arrChunks[i][loc.id][0 ..# c];
      offset += c;
    }
  }
  // NOTE: this creates a copy of states, right? :(
  var r = new BasisStates(basis.numberSpins(), states, counts, bits);
  __timer4.stop();
  return r;
}
