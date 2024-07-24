module SpaceConv

using GridapMHD
using GridapMHD: transient
using GridapPETSc, SparseMatricesCSR
using DrWatson

tf = 1.0
man_solution = [:lineartime_nonfespace,:stationary_nonfespace]
Δt = 0.5
n = 2 .^ (1:4)

params = ntuple2dict((;Δt,tf,man_solution,n))

all_params = dict_list(params)

prefix = "transient"
default = (;vtk=false)

for p in all_params
  title = savename(prefix,p)
  println("--------------------------------------")
  println("Running $title")
  println("--------------------------------------")
  transient(;title,default...,p...)
end



end # module
