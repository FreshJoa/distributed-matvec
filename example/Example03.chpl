use LatticeSymmetries;
use Time;
use RangeChunk;

proc testMultiwayMerge()
{
  var chunks: [0 .. 4, 0 .. 3] int;
  chunks[0, ..] = [45, 47,  0,  0];
  chunks[1, ..] = [0,   0,  0,  0];
  chunks[2, ..] = [4,  15, 16, 19];
  chunks[3, ..] = [17, 82,  0,  0];
  chunks[4, ..] = [0,   0,  0,  0];
  var counts = [2, 0, 4, 2, 1];
  for x in kMerge(new BlockVector(chunks, counts)) do
    writeln(x);
}

proc main() {
  initRuntime();
  defer deinitRuntime();

  testMultiwayMerge();
  var basis0 = new Basis("{\"number_spins\": 10, \"hamming_weight\": 5}");
  basis0.build();
  var buckets0 = enumerateStates(basis0, 5);
  // ref v1 = buckets0[0];
  // ref v2 = buckets0[0].toArray();
  // writeln(v1);
  // writeln(v2);
  // v2[0] = 0;
  // writeln(v1);
  // writeln(buckets0);
  const rs = determineMergeRanges(buckets0, 5);
  for (r, i) in zip(rs, 0 ..) {
    write(r, " ");
    for k in r do
      if k.size > 0 then
        write((buckets0[i, k.low], buckets0[i, k.high]), " ");
    writeln();
  }

  const blocks = statesFromHashedToBlock(buckets0);
  writeln(blocks);

  // for k in buckets0.domain {
  //   writeln(buckets0[k]);
  // }
  return 0;
}
