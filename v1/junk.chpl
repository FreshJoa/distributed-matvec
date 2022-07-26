use Search;


proc main() {
  var xs = [1, 2, 5, 9, 11, 89];
  writeln(binarySearch(xs, 2));
  writeln(binarySearch(xs, 3));
  writeln(binarySearch(xs, 10));
  writeln(binarySearch(xs, 91));
}
