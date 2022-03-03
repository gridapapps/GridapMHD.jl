
function expansion()
  title = "expansion"
  debug = false

  t = PTimer(get_part_ids(sequential,1),verbose=true)
  tic!(t,barrier=true)
  msh_file = joinpath(@__FILE__,"..","..","meshes","Expansion_coarse.msh") |> normpath
  model = GmshDiscreteModel(msh_file)
  writevtk(model,"expansion_model")
  toc!(t,"model")
  Ω = Interior(model,tags="PbLi")
  toc!(t,"triangulation")

  z = VectorValue(0,0,0)

  u_inlet((x,y,z)) = VectorValue(16*(y-1/4)*(y+1/4)*(z-1)*(z+1),0,0)

  params = Dict(
    :ptimer=>t,
    :debug=>debug,
    :fluid=>Dict(
      :domain=>Ω,
      :α=>1.0,
      :β=>1.0,
      :γ=>1.0,
      :u=>Dict(
        :tags=>["inlet", "wall"],
        :values=>[u_inlet, VectorValue(0,0,0)]
      ),
      :j=>Dict(
        :tags=>["wall"],
        :values=>[VectorValue(0,0,0)]
      ),
      :f=>z,
      :B=>z,
    ),
#    :solver=>NLSolver(show_trace=true,method=:newton,iterations=1),
  )

  xh = main(params)

  uh,ph,jh,φh = xh

  # U_u = get_fe_space(uh)
  # U_j = get_fe_space(jh)

  # tic!(t)
  # u(x) = VectorValue(sum(x),sum(x),sum(x))
  # j(x) = VectorValue(sum(x),sum(x),sum(x))
  # uh = interpolate(u,U_u)
  # jh = interpolate(j,U_j)
  # toc!(t,"interpolation")

  # eh_u = u - uh
  # eh_j = j - jh
  # dΩ = Measure(Ω,2*2)
  # eh_u_l2 = sqrt(sum(∫( eh_u⋅eh_u )dΩ))
  # eh_j_l2 = sqrt(sum(∫( eh_j⋅eh_j )dΩ))
  # toc!(t,"error_norms")

  writevtk(Ω,"$(title)_Ω",
  # cellfields=[
  # "uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,"u"=>u,"j"=>j,])
  cellfields=[
  "uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,])

  out = Dict{Symbol,Any}()
  # out[:eh_u_l2] = eh_u_l2
  # out[:eh_j_l2] = eh_j_l2

  display(t)

  out
end
