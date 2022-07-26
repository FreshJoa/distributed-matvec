module StatesIndexing {

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

}
