# sw-ising

Swendsen-Wang Algorithm for 2D square-lattice Ising model with anisotropic interaction

## Build

Requirements: CMake 3.14 or later, a C++14 compiler, Make, and Git.
The lattis dependency also requires Rust/Cargo. Install Eigen3 so CMake can find it.
An internet connection is required to download dependencies during configuration.

From the repository root, run:

```sh
cmake -S . -B build
cmake --build build
./build/ising -h
```

## Usage

Use `./build/ising -g` to print lattice information and exit without running the simulation.
The output includes site information and a list of bonds with their type, source, and target.
The periodic rectangular lattice has sizes Lx and Ly (both default to 4).
Use `-lx` and `-ly` to set them independently, for example `./build/ising -g -lx 4 -ly 6`.
Bond type 0 runs in the x direction; bond type 1 runs in the y direction.

Use `-jx Jx` and `-jy Jy` to set the interaction constants for bond types 0 (x)
and 1 (y), respectively. Both default to 1.0 and must be finite and nonnegative
(ferromagnetic interactions, with zero allowed).
The Hamiltonian is `H = -sum_b J[type(b)] s[source(b)] s[target(b)]`, with spins
`s = +/-1` and Boltzmann constant set to 1. Cluster bonds are activated with
probability `1 - exp(-2 J[type(b)] / T)` between equal spins.

For example, `./build/ising -lx 4 -ly 6 -jx 1.0 -jy 0.5 -t 2.0` runs an anisotropic model.

Licensed under the [Apache License, Version 2.0](LICENSE).
