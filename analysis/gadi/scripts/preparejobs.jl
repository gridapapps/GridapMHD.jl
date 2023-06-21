using Mustache
using DrWatson
using Dates

mkpath(datadir())

function tostr(v, n)
  lpad(v, n, '0')
end

function jobtitle(params)
  cx = params[:cx]
  px = params[:px]
  ha = params[:ha]
  nr = params[:nr]
  ns = params[:ns]
  ps = params[:ps]
  km = params[:km]
  title = "ha$(tostr(ha,5))cx$(tostr(cx,3))px$(tostr(px,2))ns$(tostr(ns,3))ps$(tostr(ps,2))km$(tostr(km,1))"
end

function petsc_options(i)
  if i == 1
    "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
  else
    error()
  end
end

function preparejobs_conv()

  # Small convergence test

  # allparams = Dict(
  #  :ha=>[500],
  #  :cx=>[15,30,45,64,90,129,179],
  #  :px=>[2],
  #  :nr=>[1],
  #  :ns=>[100],
  #  :ps=>[1],
  #  :km=>[1])
  # Take points equally spaced in a log plot 
  #log_cx=[1,1.25,1.5,1.75,2,2.25]
  #cx=map(x->Int(round(x)),10 .^ log_cx)

  allparams = Dict(
    :ha => [50],
    :cx => [10,18,32,56,100,178,400],
    :px => [2],
    :nr => [1],
    :ns => [500],
    :ps => [1],
    :km => [1])

  params = dict_list(allparams)
  dicts = map(params) do params
    cx = params[:cx]
    px = params[:px]
    ha = params[:ha]
    nr = params[:nr]
    ns = params[:ns]
    ps = params[:ps]
    km = params[:km]
    title = jobtitle(params)
    Dict(
      :q => "hugemem",
      :walltime => "04:00:00",
      :ncpus => px^2,
      :mem => "180gb",
      :jobfs => "1gb",
      :name => title,
      :n => px^2,
      :nc => (cx, cx),
      :np => (px, px),
      :B => (0.0, Float64(ha), 0.0),
      :vtk => "false",
      :title => title,
      :path => datadir(),
      :nruns => nr,
      :nsums => ns,
      :kmap => km,
      :petsc_options => petsc_options(ps),
    )
  end

  template = read(scriptsdir("jobtemplate"), String)

  for dict in dicts
    jobfile = datadir(dict[:name] * ".sh")
    open(jobfile, "w") do io
      render(io, template, dict |> tostringdict)
    end
  end

  return nothing
end

function get_walltime(tmax, tmin, n)
  msmax = Millisecond(tmax)
  msmin = Millisecond(tmin)
  ms = Millisecond(ceil(msmax.value / n))
  s = convert(Time, convert(DateTime, ms))
  smin = convert(Time, convert(DateTime, msmin))
  s < smin ? smin : s
end

# Strong Scaling
function preparejobs_ssca()
  allparams = Dict(
    :ha => [50],
    :cx => [400],
    :px => [2, 3, 4, 6, 8, 12, 16],
    :nr => [5],
    :ns => [100],
    :ps => [1],
    :km => [1])

  params = dict_list(allparams)
  dicts = map(params) do params
    cx = params[:cx]
    px = params[:px]
    ha = params[:ha]
    nr = params[:nr]
    ns = params[:ns]
    ps = params[:ps]
    km = params[:km]

    n = px^2
    nnodes = ceil(n / 48) |> Int
    title = jobtitle(params)
    ncpus = nnodes > 1 ? nnodes*48 : n

    mem=Dict(
      2=>190,
      3=>190,
      4=>190,
      6=>190,
      8=>380,
      12=>570,
      16=>950)
    que=Dict(
      2=>"normal",
      3=>"normal",
      4=>"normal",
      6=>"normal",
      8=>"normal",
      12=>"normal",
      16=>"normal")
    tw =Dict(
      2=> "08:00:00",
      3=> "07:00:00",
      4=> "06:00:00",
      6=> "05:00:00",
      8=> "04:00:00",
      12=> "03:00:00",
      16=> "02:00:00",
    )
    Dict(
      :q => que[px],
      :o => title*"o.txt",
      :e => title*"e.txt",
      :walltime => tw[px],
      :ncpus => ncpus,
      :mem => "$(mem[px])gb",
      :jobfs => "1gb",
      :name => title,
      :n => n,
      :nc => (cx, cx),
      :np => (px, px),
      :B => (0.0, Float64(ha), 0.0),
      :vtk => "false",
      :title => title,
      :path => datadir(),
      #:path => "/scratch/bt62/jp3650/data",
      :nruns => nr,
      :nsums => ns,
      :kmap => km,
      :petsc_options => petsc_options(ps),
    )
  end

  template = read(scriptsdir("jobtemplate"), String)

  for dict in dicts
    jobfile = datadir(dict[:name] * ".sh")
    open(jobfile, "w") do io
      render(io, template, dict |> tostringdict)
    end
  end

  return nothing
end


# Weak Scaling
function preparejobs_wsca()
  allparams = Dict(
    :ha => [50],
    :cx => [24],
    # :px => [2, 3, 4, 6, 8, 12, 16, 24, 32],
    :px => [2, 3, 4, 6, 8, 12, 16],
    #:px => [20],
    :nr => [6],
    :ns => [200],
    :ps => [1],
    :km => [1])

  params = dict_list(allparams)
  dicts = map(params) do params
    cx = params[:cx]
    px = params[:px]
    ha = params[:ha]
    nr = params[:nr]
    ns = params[:ns]
    ps = params[:ps]
    km = params[:km]

    n = px^2
    cx = cx * px # So N=cx^2 scales as P=px^2
    params[:cx]=cx
    nnodes = ceil(n / 48) |> Int
    title = jobtitle(params)
    Dict(
      :q => "normal",
      :o => title*"o.txt",
      :e => title*"e.txt",      
      :walltime => "04:00:00",
      :ncpus => nnodes > 1 ? nnodes*48 : n,
      :mem => "$(190*nnodes)gb",
      :jobfs => "1gb",
      :name => title,
      :n => n,
      :nc => (cx, cx),
      :np => (px, px),
      :B => (0.0, Float64(ha), 0.0),
      :vtk => "false",
      :title => title,
      :path => datadir(),
      :nruns => nr,
      :nsums => ns,
      :kmap => km,
      :petsc_options => petsc_options(ps),
    )
  end

  template = read(scriptsdir("jobtemplate"), String)

  for dict in dicts
    jobfile = datadir(dict[:name] * ".sh")
    open(jobfile, "w") do io
      render(io, template, dict |> tostringdict)
    end
  end

  return nothing
end

function preparejobs_ha_sca()
  allparams = Dict(
    :ha => [50,100,500,1000,5000,10000],
    :cx => [50,100],
    :px => [2],
    :nr => [1],
    :ns => [1000],
    :ps => [1],
    :km => [1])

  params = dict_list(allparams)
  dicts = map(params) do params
    cx = params[:cx]
    px = params[:px]
    ha = params[:ha]
    nr = params[:nr]
    ns = params[:ns]
    ps = params[:ps]
    km = params[:km]

    n = px^2
    nnodes = ceil(n / 48) |> Int
    title = jobtitle(params)
    Dict(
      :q => "normal",
      :o => title*"o.txt",
      :e => title*"e.txt",      
      :walltime => "01:00:00",
      :ncpus => nnodes > 1 ? nnodes*48 : n,
      :mem => "$(190*nnodes)gb",
      :jobfs => "1gb",
      :name => title,
      :n => n,
      :nc => (cx, cx),
      :np => (px, px),
      :B => (0.0, Float64(ha), 0.0),
      :vtk => "true",
      :title => title,
      :path => datadir(),
      :nruns => nr,
      :nsums => ns,
      :kamp_x => 1,
      :kamp_y => 1,
      :petsc_options => petsc_options(ps),
    )
  end

  template = read(scriptsdir("jobtemplate"), String)

  for dict in dicts
    jobfile = datadir(dict[:name] * ".sh")
    open(jobfile, "w") do io
      render(io, template, dict |> tostringdict)
    end
  end

  return nothing
end


