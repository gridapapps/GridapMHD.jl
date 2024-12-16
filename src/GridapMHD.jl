module GridapMHD

using Random
using LinearAlgebra
using SparseArrays, SparseMatricesCSR, BlockArrays, FillArrays
using ForwardDiff

using FileIO
using BSON
using DrWatson

using Gridap
using Gridap.Helpers, Gridap.Algebra, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField, Gridap.ODEs
using Gridap.Fields

using PartitionedArrays
using PartitionedArrays: getany

using GridapDistributed
using GridapDistributed: i_am_main

using GridapGmsh
using GridapPETSc

using GridapSolvers
using GridapSolvers.SolverInterfaces, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers
using GridapSolvers.LinearSolvers, GridapSolvers.NonlinearSolvers, GridapSolvers.BlockSolvers

# Mesh generation
include("Meshers/meshers.jl")

# Solvers
include("Solvers/gridap.jl")
include("Solvers/petsc.jl")
include("Solvers/gmg.jl")
include("Solvers/li2019.jl")
include("Solvers/badia2024.jl")

# Main driver
include("gridap_extras.jl")
include("utils.jl")
include("geometry.jl")
include("parameters.jl")
include("weakforms.jl")
include("main.jl")

# Applications
include("Applications/hunt.jl")
include("Applications/expansion.jl")
include("Applications/cavity.jl")
include("Applications/transient.jl")
include("Applications/channel.jl")
include("Applications/FullyDevelopedFlow.jl")
include("Applications/Tube.jl")
include("Applications/toy.jl")

export hunt, expansion, cavity, tube, FullyDeveloped, toy

end # module
