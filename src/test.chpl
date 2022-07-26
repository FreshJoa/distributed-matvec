use BlockDist;

proc main() {
  const Space = {1..4, 1..8};
  const targetLocales = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
  const D: domain(2) dmapped Block(boundingBox=Space, targetLocales=targetLocales) = Space;
  var A: [D] int;

  forall a in A do
    a = a.locale.id;

  writeln(A);
  return 0;
}
