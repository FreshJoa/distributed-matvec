# distributed-matvec

This repository contains Chapel implementation of a few key algorithms for
lattice-symmetries package.


## Obtaining the dependencies

[`lattice-symmetries-haskell`](https://github.com/twesterhout/lattice-symmetries-haskell)
is a required dependency. You can either compile it from source (not recommended
unless you plan to contribute patches) or download a pre-built version:

```bash
LS_HS_VERSION=abe34c4
# On Linux
wget -q https://github.com/twesterhout/lattice-symmetries-haskell/releases/download/continuous/lattice-symmetries-haskell-${LS_HS_VERSION}-ubuntu-18.04.tar.bz2
tar -xf lattice-symmetries-haskell-${LS_HS_VERSION}-ubuntu-18.04.tar.bz2
# On OS X
# wget -q https://github.com/twesterhout/lattice-symmetries-haskell/releases/download/continuous/lattice-symmetries-haskell-${LS_HS_VERSION}-macos-latest.tar.bz2
# tar -xf lattice-symmetries-haskell-${LS_HS_VERSION}-macos-latest.tar.bz2
mv lattice-symmetries-haskell-${LS_HS_VERSION} third_party
```

See also the [Github Actions workflow](.github/workflows/ci.yml) for up-to-date
installation instructions (for instance, the latest tested `LS_HS_VERSION`).


## Building

Use the [Makefile](Makefile) in this repository to build the code:

```bash
make
```

This assumes that you haven an existing installation of
[Chapel](https://chapel-lang.org/).


## Testing

To run tests, type

```bash
make check
```

This will automatically download some HDF5 input files which are required to run
the tests.


## Benchmarking

There are some benchmarks which you can run using

```bash
make benchmark-states-enumeration
make benchmark-matrix-vector-product
```

Timings should be compared with the plots from the
[`lattice-symmetries`](https://github.com/twesterhout/lattice-symmetries#bicyclist-performance)
repo.

Please, note, that running benchmarks requires some extra HDF5 input files which
will use 1-2GB of storage.


## Remarks on matrix-vector product

We want to compute `y ← H⋅x`, where `y` and `x` are distributed vectors and `H`
is an operator. `H` is never stored and is given to us in a matrix-free form.

> **A remark on notation:** in quantum mechanics it is common to use so called
> Dirac notation for vectors. So we would write |ψ⟩ when we mean that ψ is a
> vector in some (high-dimensional) vector space. In the following, I will use
> |σ⟩ do denote basis vectors of our space to avoid confusion between basis
> vectors and indices, i.e. `i` is an index (internally, an `int`) and `|σᵢ⟩` is
> the `i`'th basis vector (internally, an `uint(64)`).

Having `H` in matrix-free form means that we can efficiently compute `H|σᵢ⟩ = ∑ⱼ
cⱼ|σⱼ⟩` where `|σᵢ⟩` is a basis vector (internally represented by a bitstring,
i.e. `uint(64)`). Coefficients `cⱼ` are complex numbers `complex(128)`.
What's important to note is that we only have access to `|σⱼ⟩` rather than `j`
itself. The index `j` needs to be explicitly computed.

A simple implementation of our matrix-vector product could look like this:

```chapel
proc matrixVectorProduct(H, const ref x, ref y, const ref representatives) {
  y = 0; // initialize y with zeros
  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }
}
```

The non-trivial part is how to distribute `x`, `y`, and `representatives` among
locales. The simplest approach would be to use block-distributed vectors.
However, there are multiple drawbacks with this approach:

  * Non-uniform work distribution. The number of non-zero elements per row of
  `H` varies significantly, and if some locale only has access to the first few
  elements of `representatives` it's highly that it will have much less work to
  do than other locales.
  * ~Figuring out on which locale a particular element `|σⱼ⟩` resides, scales as
  `O(log numLocales)` which is non-optimal especially because the number of
  `|σⱼ⟩` that we have to process, is extremely big (think, 10¹⁰-10¹¹).~ This was
  wrong, the operation can actually be done in constant time if the size of
  blocks is the same across locales.

> **Note:** I'm now thinking that having a block-cyclic distribution might do
> the trick... However, the below discussion is easily transferable to the case
> when block-cyclic distribution is used instead of hashing.

We will take an alternative approach and use a simple hash function to map
states `|σ⟩` to locales:

```chapel
proc localeOf(basisState : uint(64)) : int {
  return (hash(basisState) % numLocales:uint):int;
}

// An example hash function
inline proc hash(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}
```

Unfortunately, such a distribution of states means that, unless we want to write
a custom distribution class, we have to manually manage all accesses to parts of
`x`, `y`, and `representatives` which reside on other locales.

Our `indexOf` function now looks more or less like this:

```chapel
proc indexOf(basisState : uint(64)) : int {
  const targetLocale = localeOf(basisState);
  var index : int;
  on Locales[targetLocale] {
    index = localIndexOf(basisState);
  }
  return index;
}
```

This is, of course, terribly inefficient because for each `|σⱼ⟩` we make a
remote procedure call (which cannot be optimized because of the `localIndexOf`
call).

Instead, we can rewrite the inner loop of our algorithm as follows:

```chapel
for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
  const t = x[i] * cⱼ;
  const targetLocale = localeOf(|σⱼ⟩);
  on targetLocale {
    // copy stuff to targetLocale
    const t' = t;
    const |σⱼ⟩' = |σⱼ⟩;
    // everything else is done locally
    const j = localIndexOf(|σⱼ⟩');
    y.localAccess(j) += t';
  }
}
```

This is an improvement, because we now go from `here` to `targetLocale`, but
never have to communicate anything back. However, we still have the problem that
a new task is spawned for each `|σⱼ⟩`. Obvious solution is to batch together
many `|σⱼ⟩`s and then process them all at once:

```chapel
for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
  const t = x[i] * cⱼ;
  const targetLocale = localeOf(|σⱼ⟩);
  magic.enqueue(targetLocale, |σⱼ⟩, t);
}
```

Here, `magic` record has a buffer for each `targetLocale` which it fills up.
When the buffer is full, `magic` will spawn a task on `targetLocale` to process
it:

```chapel
proc process(targetLocale : locale, buffer : [] (uint(64), complex(128))) {
  on targetLocale {
    const localBuffer = buffer;
    for (t, |σⱼ⟩) in localBuffer {
      const j = localIndexOf(|σⱼ⟩);
      y.localAccess(j) += t;
    }
  }
}
```

There is a problem with the above approach, though, which is that we have a race
condition: multiple tasks might try to perform `y.localAccess(j) += t` in
parallel. If `y.eltType` was an integral type, we could have used atomic
operations, but since `y.eltType` is `real` or `complex`, we need some kind of
locking mechanism. The simplest is to protect `y` with a mutex or a spinlock,
but that does not scale to many cores. Instead, we can have an array of locks
protecting parts of `y`. Then, since accesses to `y` are pretty random, there is
hope that different tasks will access different parts of `y` and won't have to
wait:

```chapel
record ConcurrentAccessor {
  type eltType;
  var _data : c_ptr(eltType);
  var _numElts : int;
  var _blockSize : int;
  var numLocks : int;
  var _locks : [0 ..# numLocks] chpl_LocalSpinlock;

  proc init(ref arr : [] ?t, numLocks : int = concurrentAccessorNumLocks)
      where !arr.domain.stridable && arr.domain.rank == 1 {
    this.eltTYpe = t;
    this._data = c_ptrTo(arr[arr.domain.low]);
    this._numElts = arr.size;
    this._blockSize = (arr.size + numLocks - 1) / numLocks;
    this.numLocks = numLocks;
    assert(_blockSize * numLocks >= _numElts);
  }

  inline proc _blockIdx(i : int) : int {
    return i / _blockSize;
  }

  inline proc localAdd(i : int, x : eltType) {
    const blockIdx = _blockIdx(i);
    ref lock = _locks[blockIdx];
    lock.lock();
    _data[i] += x;
    lock.unlock();
  }
}
```

Now, instead of doing `y.localAccess(i) += t` we can do

```
const yAccessor = new ConcurrentAccessor(y);
...
yAccessor.localAdd(i, t); // similar to y[i] += t, but thread-safe
```












