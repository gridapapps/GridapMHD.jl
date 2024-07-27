module TransientConv

using GridapMHD
using GridapMHD: transient
using GridapPETSc, SparseMatricesCSR
using DrWatson

n = 2
tf = 1.0
man_solution = [:nonlineartime_fespace]
Δt = 2. .^ -(0:4)

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
