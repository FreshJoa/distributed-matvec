module MultiwayMerge {

use FFI;
use Vector;

use List;
use Search;
use CTypes;
use DynamicIters;

record Node {
  var chunkIndex : int;
  var indexInChunk : int;
  var outerParentIndex : int;

  inline proc isEmpty() { return chunkIndex == -1; }
  inline proc isInfinity() { return indexInChunk == -1; }
}

private inline proc indexAccess(const ref chunks, i : int, j) const ref {
  return chunks[i][j];
}
private inline proc indexAccess(const ref chunks : BlockVector, i : int, j) const ref {
  return chunks[i, j];
}
private inline proc indexAccess(const ref chunks : [?D] ?eltType, i : int, j) const ref
    where D.rank == 2 {
  return chunks[i, j];
}
private inline proc indexAccess(const ref chunks : [] Vector(?t), i : int, j) const ref
    where chunks.domain.rank == 1 {
  return chunks[i][j];
}
private inline proc indexAccess(const ref chunks : [] (c_ptr(?t), int), i : int, j) const ref
    where chunks.domain.rank == 1 {
  const (ptr, size) = chunks[i];
  assert(j < size);
  return ptr[j];
}

private inline proc chunkSize(const ref chunks : BlockVector, i : int) {
  return chunks.count(i);
}
private inline proc chunkSize(const ref chunks : [] Vector(?t), i : int)
    where chunks.domain.rank == 1 {
  return chunks[i].size;
}
private inline proc chunkSize(const ref chunks : [] ?t, i : int)
    where chunks.domain.rank == 2 {
  return chunks.dim(1).size;
}
private inline proc chunkSize(const ref chunks : [] (c_ptr(?t), int), i : int)
    where chunks.domain.rank == 1 {
  const (ptr, size) = chunks[i];
  return size;
}

private inline proc numChunks(const ref chunks : BlockVector) {
  return chunks.numBlocks;
}
private inline proc numChunks(const ref chunks : [] Vector(?t))
    where chunks.domain.rank == 1 {
  return chunks.size;
}
private inline proc numChunks(const ref chunks : [] ?t)
    where chunks.domain.rank == 2 {
  return chunks.dim(0).size;
}
private inline proc numChunks(const ref chunks : [] (c_ptr(?t), int))
    where chunks.domain.rank == 1 {
  return chunks.size;
}

private inline proc _isLess(const ref a : Node, const ref b : Node,
                            const ref chunks) : bool {
  assert(!a.isEmpty() && !b.isEmpty());
  if a.isInfinity() { return false; }
  if b.isInfinity() { return true; }
  return indexAccess(chunks, a.chunkIndex, a.indexInChunk)
           < indexAccess(chunks, b.chunkIndex, b.indexInChunk);
}

private inline proc _makeNode(chunkIndex : int, in indexInChunk : int,
                              outerParentIndex : int, chunks) {
  if chunkIndex >= numChunks(chunks) || indexInChunk >= chunkSize(chunks, chunkIndex) then
    indexInChunk = -1;
  return new Node(chunkIndex, indexInChunk, outerParentIndex);
}


record TournamentTree {
  var numberChunks : int;
  var nodes : [0 ..# numberChunks] Node;
  var chunks;

  proc init(chunks)
  {
    this.numberChunks = numChunks(chunks) + (numChunks(chunks) % 2);
    this.nodes = [i in 0 ..# numberChunks] new Node(-1, -1, -1);
    this.chunks = c_const_ptrTo(chunks);
  }

  inline proc innerParentIndex(i : int) {
    assert(i > 0);
    return i / 2;
  }

  proc insert(i : int, const ref node : Node) : void {
    assert(i != -1);
    assert(!node.isEmpty());
    if nodes[i].isEmpty() { nodes[i] = node; }
    else {
      const parent : int = innerParentIndex(i);
      if _isLess(node, nodes[i], chunks.deref()) {
        insert(parent, node);
      }
      else {
        const old = nodes[i];
        nodes[i] = node;
        insert(parent, old);
      }
    }
  }
  inline proc insert(const ref node : Node) : void {
    insert(node.outerParentIndex, node);
  }

  inline proc processOne() {
    const result = nodes[0];
    if (result.isInfinity()) { return result; }
    nodes[0] = new Node(-1, -1, -1);
    insert(_makeNode(result.chunkIndex, result.indexInChunk + 1,
                     result.outerParentIndex, chunks.deref()));
    return result;
  }
}

private proc _buildTree(const ref chunks) {
  var tree = new TournamentTree(chunks);
  const numberChunks = tree.numberChunks;
  for i in 0 .. numberChunks / 2 - 1 {
    const outerParentIndex = numberChunks - 1 - i;
    tree.insert(_makeNode(numberChunks - 1 - 2 * i, 0, outerParentIndex, chunks));
    tree.insert(_makeNode(numberChunks - 2 - 2 * i, 0, outerParentIndex, chunks));
  }
  return tree;
}

iter kMergeIndices(chunks) {
  var tree = _buildTree(chunks);
  var root = tree.processOne();
  while (!root.isInfinity()) {
    yield (root.chunkIndex, root.indexInChunk);
    root = tree.processOne();
  }
}

iter kMerge(chunks) {
  for (chunkIndex, indexInChunk) in kMergeIndices(chunks) {
    yield indexAccess(chunks, chunkIndex, indexInChunk);
  }
}

/*
var _chunkBoundsTime = new MeasurementTable("Merge.chunkBounds");

private proc chunkBounds(const ref arrays, const ref counts, chunkSize : int) {
  assert(chunkSize >= 1);
  assert(&& reduce [i in counts.domain] counts[i] > 0);
  var __time = getTimerFor(_chunkBoundsTime);
  type eltType = indexAccess(arrays, 0, 0).type;
  const lower = indexAccess(arrays, 0, 0);
  const upper = indexAccess(arrays, 0, counts[0] - 1);
  // const lower = min reduce [i in counts.domain] indexAccess(arrays, i, 0);
  // const upper = max reduce [i in counts.domain] indexAccess(arrays, i, counts[i] - 1);

  var bounds : list(eltType);
  bounds.append(lower);
  var offset = chunkSize;
  while (offset < counts[0]) {
    bounds.append(indexAccess(arrays, 0, offset));
    offset += chunkSize;
  }
  if (bounds.last() != upper || bounds.size < 2) { bounds.append(upper); }
  // writeln("chunkBounds: ", bounds);
  var boundsArr : [0 ..# bounds.size] eltType = bounds;
  return boundsArr;
}

var _chunkOffsetsTime = new MeasurementTable("Merge.chunkOffsets");

private proc chunkOffsets(const ref arrays, const ref counts, chunkSize : int) {
  var __time = getTimerFor(_chunkOffsetsTime);
  const bounds = chunkBounds(arrays, counts, chunkSize);

  var offsets : [0 ..# counts.size, 0 ..# bounds.size] int;
  // startVerboseComm();
  forall i in counts.domain {
    const localBounds = bounds;
    const localCount = counts[i];
    var localOffsets : [0 ..# localBounds.size] int;
    localOffsets[0] = 0;
    localOffsets[localBounds.size - 1] = localCount;
    for j in 1 .. bounds.size - 2 {
      if arrays.domain.rank == 1 {
        localOffsets[j] = binarySearch(
          arrays[i][localOffsets[j - 1] .. counts[i] - 1], localBounds[j])[1];
      }
      else {
        localOffsets[j] = binarySearch(
          arrays[i, localOffsets[j - 1] .. counts[i] - 1], localBounds[j])[1];
      }
    }
    offsets[i, ..] = localOffsets;
  }
  // stopVerboseComm();
  return offsets;
}

private inline proc linearizeOffsets(const ref offsets) {
  return [j in offsets.dim(1)] (+ reduce offsets[.., j]);
}

var _chunkedIterBodyTime = new MeasurementTable("Merge.chunkedIterBody");
var _kMergeIndicesChunked = new MeasurementTable("Merge.kMergeIndicesChunked");

private proc chunkedIterBody(chunkId : int, const ref arrays, const ref counts,
                             const ref offsets, numberArrays, type eltType) {
  var __time = getTimerFor(_chunkedIterBodyTime);
  // Gather chunk
  const localCounts : [0 ..# numberArrays] int =
    [i in 0 ..# numberArrays] offsets[i, chunkId + 1] - offsets[i, chunkId];
  const totalLocalCount = + reduce localCounts;
  const maxLocalCount = max reduce localCounts;
  var localArrays : [0 ..# numberArrays, 0 ..# maxLocalCount] eltType;
  foreach i in 0 ..# numberArrays {
    if arrays.domain.rank == 1 {
      localArrays[i, 0 ..# localCounts[i]] =
        arrays[i][offsets[i, chunkId] .. offsets[i, chunkId + 1] - 1];
    }
    else {
      localArrays[i, 0 ..# localCounts[i]] =
        arrays[i, offsets[i, chunkId] .. offsets[i, chunkId + 1] - 1];
    }
  }
  // writeln("chunkId=", chunkId, " localArrays:");
  // writeln(localArrays);
  var localMergeIndices : [0 ..# totalLocalCount] (int, int);
  for ((arrayIndex, indexInArray), j) in zip(kMergeIndices(localArrays, localCounts),
                                             localMergeIndices) {
    const o = offsets[arrayIndex, chunkId];
    j = (arrayIndex, o + indexInArray);
  }
  return localMergeIndices;
}

iter kMergeIndicesChunked(const ref arrays, const ref counts, chunkSize : int) {
  var __time = getTimerFor(_kMergeIndicesChunked);
  const offsets = chunkOffsets(arrays, counts, chunkSize);
  const linear = linearizeOffsets(offsets);
  const numberArrays = counts.size;
  const numberChunks = offsets.dim(1).size - 1;
  type eltType = indexAccess(arrays, 0, 0).type;
  for chunkId in 0 ..# numberChunks {
    const indices = chunkedIterBody(chunkId, arrays, counts, offsets, numberArrays, eltType);
    yield (linear[chunkId], indices);
  }
}
iter kMergeIndicesChunked(param tag : iterKind,
                          const ref arrays, const ref counts, chunkSize : int)
    where tag == iterKind.standalone {
  var __time = getTimerFor(_kMergeIndicesChunked);
  const offsets = chunkOffsets(arrays, counts, chunkSize);
  // const numberArrays = counts.size;
  // const numberChunks = offsets.dim(1).size - 1;
  type eltType = indexAccess(arrays, 0, 0).type;
  // writeln("[Chapel] using parallel version of kMergeIndicesChunked; ",
  //         numberChunks, " to process...");
  coforall loc in Locales do on loc {
    const localOffsets = offsets;
    const linear = linearizeOffsets(localOffsets);
    const numberArrays = counts.size;
    const numberChunks = localOffsets.dim(1).size - 1;
    const chunkIds = 0 ..# numberChunks by numLocales align loc.id;
    forall chunkId in chunkIds {
      const indices = chunkedIterBody(chunkId, arrays, counts,
                                      localOffsets, numberArrays, eltType);
      yield (linear[chunkId], indices);
    }
  }
  // forall chunkId in 0 ..# numberChunks {
  //   const indices = chunkedIterBody(chunkId, arrays, counts, offsets, numberArrays, eltType);
  //   yield (linear[chunkId], indices);
  // }
  writeln("[Chapel] Done with kMergeIndicesChunked!");
}
*/

} // end module
