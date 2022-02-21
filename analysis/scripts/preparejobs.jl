using Mustache
using DrWatson
using Dates

mkpath(datadir())

function tostr(v,n)
  lpad(v,n,'0')
end

function jobtitle(params)
  cx = params[:cx]
  px = params[:px]
  ha = params[:ha]
  nr = params[:nr]
  ns = params[:ns]
  ps = params[:ps]
  title="ha$(tostr(ha,5))cx$(tostr(cx,3))px$(tostr(px,2))ns$(tostr(ns,3))ps$(tostr(ps,2))"
end

function petsc_options(i)
  if i == 1
    "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
  else
    error()
  end
end

allparams = Dict(
 :ha=>[500],
 :cx=>[15,30,45,64,90,129,179],
 :px=>[2],
 :nr=>[1],
 :ns=>[100],
 :ps=>[1])

params = dict_list(allparams)
dicts = map(params) do params
  cx = params[:cx]
  px = params[:px]
  ha = params[:ha]
  nr = params[:nr]
  ns = params[:ns]
  ps = params[:ps]
  title=jobtitle(params)
  Dict(
   :q=>"normal",
   :walltime=>"00:30:00",
   :ncpus=>px^2,
   :mem=>"180gb",
   :jobfs=>"2gb",
   :name=>title,
   :n=>px^2,
   :nc=>(cx,cx),
   :np=>(px,px),
   :B=>(0.,Float64(ha),0.),
   :debug=>"false",
   :vtk=>"false",
   :title=>title,
   :path=>datadir(),
   :nruns=>nr,
   :nsums=>ns,
   :petsc_options=>petsc_options(ps),
  )
end

template = read(scriptsdir("jobtemplate"),String)

for dict in dicts
  jobfile = datadir(dict[:name]*".sh")
  open(jobfile,"w") do io
    render(io,template,dict|>tostringdict)
  end
end


