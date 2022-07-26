module BasisLookupCache {


record LookupCache {
  type eltType = uint(64);
  var _count : int;
  var _representatives : [0 ..# _count] eltType;
}

}
