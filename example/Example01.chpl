use LatticeSymmetries;
use Time;

proc main() {
  initRuntime();
  defer deinitRuntime();

  var basis = SpinBasis(4, 2);
  writeln(basis.numberSites());
  writeln(basis.isHammingWeightFixed());

  basis = SpinBasis("{\"number_spins\": 10, \"hamming_weight\": 5, \"spin_inversion\": -1,\
                      \"symmetries\": [\
                        {\"permutation\": [1, 2, 3, 4, 5, 6, 7, 8, 9, 0], \"sector\": 5},\
                        {\"permutation\": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0], \"sector\": 1}\
                      ] }");
  var tuples = reshape(
      [0, 1,
       1, 2,
       2, 3,
       3, 4,
       4, 5,
       5, 6,
       6, 7,
       7, 8,
       8, 9,
       9, 0], {0 ..# 10, 0 ..# 2});
  var a = new Operator(basis, "σᶻ₀ σᶻ₁", tuples); 
  var b = new Operator(basis, "σˣ₀ σˣ₁", tuples);
  var c = new Operator(basis, "σʸ₀ σʸ₁", tuples);
  var timer = new Timer();
  timer.start();
  var op = a + b + c;
  timer.stop();
  writeln(timer.elapsed());

  timer.clear();
  timer.start();
  basis.build();
  timer.stop();
  writeln("building representatives took ", timer.elapsed());
  // writeln(localEnumerateRepresentatives(basis, 0, 16));

  ref states = basis.representatives();
  // var x : [0 ..# 16] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966, 0.24723682304851302, 0.38044071869638985, 0.053766497354555076, 0.1469359882375053, 0.5787416396568855, 0.951776938627626, 0.5199039424498677, 0.7953057411693661, 0.9261507114991209, 0.700102655754242];
  // var y : [0 ..# 16] real(64);
  // var x : [0 ..# 6] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966];
  // var y : [0 ..# 6] real(64);
  var x : [0 ..# 16] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966, 0.24723682304851302, 0.38044071869638985, 0.053766497354555076, 0.1469359882375053, 0.5787416396568855, 0.951776938627626, 0.5199039424498677, 0.7953057411693661, 0.9261507114991209, 0.700102655754242];
  var y : [0 ..# 16] real(64);

  timer.clear();
  timer.start();
  localMatVec(op, states, x, y);
  timer.stop();
  writeln(y);
  writeln(timer.elapsed());
}
