import h5py
import lattice_symmetries as ls
import numpy as np
import os
import time
import yaml

np.random.seed(42)


def load_hamiltonian(filename: str):
    _, extension = os.path.splitext(filename)
    if extension != ".yaml" and extension != ".yml":
        raise ValueError(
            "Do not know how to read the basis from files with extension '{}'..."
            " Try specifying a path to a YAML file instead.".format(extension)
        )

    with open(filename, "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    basis = ls.SpinBasis.load_from_yaml(config["basis"])
    return ls.Operator.load_from_yaml(config["hamiltonian"], basis)


def generate(basis_filename: str, output_filename: str, batch_size: int = 1, save: bool = True):
    hamiltonian = load_hamiltonian(basis_filename)
    hamiltonian.basis.build()
    x = np.random.rand(hamiltonian.basis.number_states, batch_size) - 0.5
    x = np.asfortranarray(x)
    tick = time.time()
    y = hamiltonian(x)
    tock = time.time()
    print("{} spent in matrix-vector using OpenMP".format(tock - tick))

    if save:
        with h5py.File(output_filename, "w") as out:
            out["/x"] = x.T
            out["/y"] = y.T


def main():
    # On cn18 (32-core Intel Xeon):
    #   - 24 sites:  0.0928194522857666 spent in matrix-vector using OpenMP
    #   - 28 sites:  1.5801198482513428 spent in matrix-vector using OpenMP
    #   - 30 sites:  7.331966161727905 spent in matrix-vector using OpenMP
    #   - 32 sites:  34.43167591094971 spent in matrix-vector using OpenMP
    # On cn20 (20-core Intel Xeon):
    #   - 24 sites: 0.17744112014770508 spent in matrix-vector using OpenMP
    #   - 28 sites: 2.703888177871704 spent in matrix-vector using OpenMP
    #   - 30 sites: 12.458613872528076 spent in matrix-vector using OpenMP
    #   - 32 sites: 55.81617712974548 spent in matrix-vector using OpenMP
    # On cn71 (64-core AMD Epyc):
    #   - 24 sites: 0.040898799896240234 spent in matrix-vector using OpenMP
    #   - 28 sites: 0.6611583232879639 spent in matrix-vector using OpenMP
    #   - 30 sites: 2.513392210006714 spent in matrix-vector using OpenMP
    #   - 32 sites: 11.471940994262695 spent in matrix-vector using OpenMP
    for i in [4, 6, 8, 10, 12, 16, 20, 24, 28, 30, 32]:
        prefix = "large-scale/" if i >= 24 else ""
        generate(
            "data/old/heisenberg_chain_{}.yaml".format(i),
            "data/{}matvec/heisenberg_chain_{}.h5".format(prefix, i),
            1,
            save=True,
        )


if __name__ == "__main__":
    main()
