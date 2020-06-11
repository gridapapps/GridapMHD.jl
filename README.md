# GridapMHD
[![Run CI](https://github.com/gridapapps/GridapMHD.jl/workflows/Run%20CI/badge.svg?branch=master)](https://github.com/gridapapps/GridapMHD.jl/actions?query=workflow%3A%22Run+CI%22)
[![codecov](https://codecov.io/gh/gridapapps/GridapMHD.jl/branch/master/graph/badge.svg?token=eSKW5MrXYz)](https://codecov.io/gh/gridapapps/GridapMHD.jl)

This application repository includes an incompressible inductionless MHD FE
formulation. It also include a series of tools to facilitate building new
[Gridap](https://github.com/gridap/Gridap.jl) drivers for magnetohydrodynamics.



## Installation
Enter Julia package manager (type `]`). Install the master branch of this repository as follows:
```
pkg> add https://github.com/gridapapps/GridapMHD.jl.git
```

## Usage
Include the drivers running
```
julia> using GridapMHD:<drivername>
```
Then execute them by calling the main function with the following syntax:
```
julia> <drivername>(optional_driver_kwargs...)
```
This function returns a tuple containing the solution as a `Gridap.MultiFieldFEFunction`,
the triangulation used, and the quadrature.

Current implemented drivers are:
- Shercliff
- Hunt

## Minimal Example
Minimal working example:
```
pkg> add https://github.com/gridapapps/GridapMHD.jl
julia> using GridapMHD:hunt
julia> xh, trian, quad = hunt(nx=3,ny=3,resultsfile="results.vtu");
julia> uh, ph, jh, φh = xh
```
This will solve Hunt's problem in Ω=[-1,1]x[-1,1]x[0,0.1] with a 3x3x3 mesh.
(see [this](https://www.sciencedirect.com/science/article/pii/S0021999111000088) for details of Hunt's problem setting). By default a Re=10 and
Ha=10 are used. By default no result file is generated, but if `resultsfile=
<filename>` is specified, the solution will be saved in `<filename>` VTK file.

## How to cite
If you have used these drivers in a scientific publication, please cite Gridap library as follows:

```
@article{gridap_guide_2019,
    author={Francesc Verdugo and Santiago Badia},
    journal = {{arXiv}},
    title = {{A user-guide to Gridap -- grid-based approximation of partial differential equations in Julia}},
    year = {2019},
    eprint={1910.01412},
    archivePrefix={arXiv},
    primaryClass={cs.MS},
}
```
