# GridapMHD
[![Run CI](https://github.com/gridapapps/GridapMHD.jl/workflows/Run%20CI/badge.svg?branch=master)](https://github.com/gridapapps/GridapMHD.jl/actions?query=workflow%3A%22Run+CI%22)
[![codecov](https://codecov.io/gh/gridapapps/GridapMHD.jl/branch/master/graph/badge.svg?token=eSKW5MrXYz)](https://codecov.io/gh/gridapapps/GridapMHD.jl)

This package implements a finite element solver for the steady-state, incompressible, and inductionless magneto-hydro-dinamics (MHD) problem.
The formulation is taken from [DOI:10.1137/19M1260372](https://doi.org/10.1137/19M1260372) and it is implemented using the tools provided by [Gridap](https://github.com/gridap/Gridap.jl), a free and open-source finite element library fully implemented in the Julia programming language.

## Installation
Enter Julia package manager (type `]`). Install the master branch of this repository as follows:
```
pkg> add https://github.com/gridapapps/GridapMHD.jl.git
```

## Usage

The library provides a driver function `GridapMHD.main` that takes a dictionary containing several parameters defining a MHD problem and returns an object representing the solution of the MHD problem.

```julia
using GridapMHD
params = Dict(...)
u,p,j,φ = GridapMHD.main(params)
```

The returned value is of type `Gridap.MultiField.MultiFieldFEFunction`, which can be unpacked to get access to the different fields of the MHD solution (fluid velocity, fluid pressure, charge current, and electric potential respectively). One can further post process these quantities using the tools provided by Gridap.

## Minimal Example

A demo example is provided in function `GridapMHD.hunt`, which implements a Hunt benchmark. This function prepares the `params` dictionary for this particular case and calls `GridapMHD.main(params)` to solve the problem. Consider function `GridapMHD.hunt` as a starting point for preparing more complex computations.

## How to cite
If you have used these drivers in a scientific publication, please cite Gridap library as follows:

```
@article{Badia2020,
  doi = {10.21105/joss.02520},
  url = {https://doi.org/10.21105/joss.02520},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2520},
  author = {Santiago Badia and Francesc Verdugo},
  title = {Gridap: An extensible Finite Element toolbox in Julia},
  journal = {Journal of Open Source Software}
}
```
