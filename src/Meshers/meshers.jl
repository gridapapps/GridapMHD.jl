module Meshers
  using DrWatson
  using PartitionedArrays
  using Gridap, GridapP4est, GridapSolvers

  using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.Adaptivity

  include("expansion_mesher.jl")

  export expansion_generate_base_mesh
  export expansion_generate_mesh
  export expansion_generate_mesh_hierarchy
  
end