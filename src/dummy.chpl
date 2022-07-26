use RangeChunk;
use CTypes;
use Time;
use Random;

private inline proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

inline proc localeIdxOf(basisState : uint(64)) : int {
  return (hash64_01(basisState) % numLocales:uint):int;
}

record PartitionInfo {
    var _countOrOffset : int;
    var nextOffset : int;

    inline proc count ref { return _countOrOffset; }
    inline proc offset ref { return _countOrOffset; }
};

proc partitionBy(in first : c_ptr(?eltType), last : c_ptr(eltType), predicate) {
    while true {
      if first == last then return last;
      if !predicate(first.deref()) then break;
      first += 1;
    }

    var it = first + 1;
    while it != last {
      if predicate(it.deref()) {
        first.deref() <=> it.deref();
        first += 1;
      }
      it += 1;
    }
    return first;
}

proc radixOneStep(keys : [] uint(8), offsets : [] int, arrs...?numArrs)
{
  var partitions : c_array(PartitionInfo, 256);
  foreach key in keys {
    partitions[key:int].count += 1;
  }

  var remainingPartitions : c_array(uint(8), 256);
  var numPartitions : int;
  var total : int;
  for i in 0 ..# 256 {
    const count = partitions[i].count;
    if count > 0 {
      partitions[i].offset = total;
      total += count;
      remainingPartitions[numPartitions] = i:uint(8);
      numPartitions += 1;
    }
    partitions[i].nextOffset = total;
  }

  var lastRemaining = remainingPartitions:c_ptr(uint(8)) + numPartitions;
  var endPartition = remainingPartitions:c_ptr(uint(8)) + 1;
  while lastRemaining - endPartition > 0 {
    record Func {
      inline proc this(partitionIdx : uint(8)) {
        ref beginOffset = partitions[partitionIdx:int].offset;
        ref endOffset = partitions[partitionIdx:int].nextOffset;
        if beginOffset == endOffset then return false;

        for i in beginOffset .. endOffset - 1 {
          ref offset = partitions[keys[i]:int].offset;
          keys[i] <=> keys[offset];
          foreach k in 0 ..# numArrs {
            ref arr = arrs[k];
            arr[i] <=> arr[offset];
          }
          offset += 1;
        }
        return beginOffset != endOffset;
      }
    }
    lastRemaining = partitionBy(remainingPartitions:c_ptr(uint(8)), lastRemaining, new Func());
  }

  offsets[0] = 0;
  foreach i in 1 ..# 256 {
    offsets[i] = partitions[i - 1].nextOffset;
  }
}

inline proc GET(addr, node, rAddr, size) {
  __primitive("chpl_comm_get", addr, node, rAddr, size);
}

inline proc PUT(addr, node, rAddr, size) {
  __primitive("chpl_comm_put", addr, node, rAddr, size);
}

config const bufferSize = 1000000;

proc main() {
  // var arr : [0 ..# 5] uint(64) = 1:uint .. 5:uint;
  // var keys : [0 ..# 5] uint(8) = [x in arr] localeIdxOf(x):uint(8);
  // writeln(arr);
  // writeln(keys);
  // var counts : [0 ..# 256 + 1] int;
  // radixOneStep(keys, counts, arr);
  // writeln(arr);
  // writeln(keys);
  // writeln(counts[0 ..# numLocales + 1 + 5]);

  var arr : [0 ..# bufferSize] real;
  fillRandom(arr);
  const arrPtr = c_ptrTo(arr[0]);

  coforall loc in Locales do on loc {
    var myArr : [0 ..# bufferSize] real;
    var timer = new Timer();
    timer.start();
    GET(c_ptrTo(myArr[0]), 0, arrPtr, bufferSize:c_size_t * c_sizeof(real));
    timer.stop();
    writeln("Sum: ", + reduce myArr, "; time: ", timer.elapsed());
  }
}
