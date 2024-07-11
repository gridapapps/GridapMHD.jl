module GridapMHD

using Random
using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using BlockArrays
using ForwardDiff

using FileIO
using BSON
using DrWatson

using Gridap
using Gridap.Helpers, Gridap.Algebra, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField, Gridap.ODEs

using PartitionedArrays
using PartitionedArrays: getany

using GridapDistributed
using GridapDistributed: i_am_in, i_am_main

using GridapGmsh
using GridapPETSc

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.NonlinearSolvers, GridapSolvers.BlockSolvers

# Mesh generation
include("Meshers/meshers.jl")

# Solvers
include("Solvers/gridap.jl")
include("Solvers/petsc.jl")
include("Solvers/li2019.jl")

# Main driver
include("Fixes.jl")
include("ExtraFunctions.jl")
include("parameters.jl")
include("weakforms.jl")
include("Main.jl")

# Applications
include("Applications/hunt.jl")
include("Applications/expansion.jl")
include("Applications/cavity.jl")
include("Applications/transient.jl")
include("Applications/channel.jl")
include("Applications/FullyDevelopedFlow.jl")
include("Applications/Tube.jl")

export hunt, expansion, cavity, tube, FullyDeveloped

end # module
