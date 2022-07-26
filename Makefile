.RULES:
.POSIX:

UNAME = $(shell uname)
OPTIMIZATION ?= --debug
COMPILER ?= llvm
CFLAGS = -Ithird_party/include $(OPTIMIZATION) --target-compiler=$(COMPILER)
LDFLAGS += -Lthird_party/lib -llattice_symmetries_haskell -llattice_symmetries_core
# ifeq ($(UNAME), Linux)
#   LDFLAGS += -lnuma
# endif

PRIMME_CFLAGS = -I/home/tom/src/primme/include
PRIMME_LDFLAGS = -L/home/tom/src/primme/lib -lprimme -lopenblas -lm -lgomp -lpthread

# MODULES = src/ApplyOperator.chpl src/StatesEnumeration.chpl src/helper.c
MODULES = src/LatticeSymmetries.chpl \
	  src/FFI.chpl \
	  src/HDF5.chpl \
	  src/ForeignTypes.chpl \
	  src/StatesEnumeration.chpl \
	  src/ConcurrentAccessor.chpl \
	  src/BatchedOperator.chpl \
	  src/CommunicationQueue.chpl \
	  src/DistributedMatrixVector.chpl \
	  src/MultiwayMerge.chpl \
	  src/Vector.chpl \
	  src/helper.c

.PHONY: all
all: examples

ifeq ($(UNAME), Linux)
  CHPL_LIBS = LD_LIBRARY_PATH=$(PWD)/third_party/lib:$$LD_LIBRARY_PATH
endif
ifeq ($(UNAME), Darwin)
  CHPL_LIBS = DYLD_LIBRARY_PATH=$(PWD)/third_party/lib:$$DYLD_LIBRARY_PATH
endif
CHPL_ARGS =

.PHONY: examples
examples: bin/Example02 bin/Example05

.PHONY: test
test: bin/TestStatesEnumeration bin/TestMatrixVectorProduct

.PHONY: check
check: check-states-enumeration check-matrix-vector-product

.PHONY: benchmark-states-enumeration
benchmark-states-enumeration: bin/TestStatesEnumeration data/large-scale
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/large-scale/construction/heisenberg_pyrochlore_2x2x2.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/large-scale/construction/heisenberg_kagome_36.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_square_6x6.yaml --kRepresentatives data/large-scale/construction/heisenberg_square_6x6.h5

.PHONY: check-states-enumeration
check-states-enumeration: bin/TestStatesEnumeration data/construction
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_chain_10.yaml --kRepresentatives data/construction/heisenberg_chain_10.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_12.yaml --kRepresentatives data/construction/heisenberg_kagome_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_16.yaml --kRepresentatives data/construction/heisenberg_kagome_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/old/heisenberg_square_4x4.yaml --kRepresentatives data/construction/heisenberg_square_4x4.h5

.PHONY: benchmark-matrix-vector-product
benchmark-matrix-vector-product: bin/TestMatrixVectorProduct data/large-scale
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/heisenberg_pyrochlore_2x2x2.h5
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/heisenberg_kagome_36.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_6x6.yaml --kVectors data/large-scale/matvec/heisenberg_square_6x6.h5

.PHONY: check-matrix-vector-product
check-matrix-vector-product: bin/TestMatrixVectorProduct data/matvec
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_4.yaml --kVectors data/matvec/heisenberg_chain_4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_10.yaml --kVectors data/matvec/heisenberg_chain_10.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_16.yaml --kVectors data/matvec/heisenberg_kagome_16.h5


TEST_DATA_URL = https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download

data/construction:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Fconstruction && \
	unzip tmp.zip && rm tmp.zip

.PHONY: data/matvec
data/matvec:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Fmatvec && \
	unzip tmp.zip && rm tmp.zip

.PHONY: data/large-scale
data/large-scale:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Flarge-scale && \
	unzip tmp.zip && rm tmp.zip

lib: lib/liblattice_symmetries_chapel.so

lib/liblattice_symmetries_chapel.so: $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) --library --dynamic -o lattice_symmetries_chapel $^ $(LDFLAGS)

bin/TestStatesEnumeration: test/TestStatesEnumeration.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/TestMatrixVectorProduct: test/TestMatrixVectorProduct.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/Example01: example/Example01.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/Example02: example/Example02.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/Example03: example/Example03.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/Example04: example/Example04.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/Example05: example/Example05.chpl $(MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ --main-module $(@F) $^ $(LDFLAGS)

bin/dummy: src/dummy.chpl
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ $<


# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
.PHONY: error
error: src/error.chpl
	chpl -o $@ $^

.PHONY: clean
clean:
	rm -rf bin/
