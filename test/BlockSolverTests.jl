
using Gridap, GridapDistributed, GridapSolvers
using Gridap.MultiField
using Gridap.Algebra

using GridapMHD: _hunt, add_default_params, _fluid_mesh, weak_form, _find_backend, p_conformity

function hunt(;
  backend=nothing,
  np=nothing,
  parts=nothing,
  title = "hunt",
  path=".",
  kwargs...)

  @assert parts === nothing
  if backend === nothing
    @assert np === nothing
    return _hunt(;title=title,path=path,kwargs...)
  else
    @assert backend !== nothing
    return with_backend(_find_backend(backend),(np...,1)) do _parts
      _hunt(;parts=_parts,title=_title,path=path,kwargs...)
    end
  end
end

_params = hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=false,
  solver=:block_cg,
)

params = add_default_params(_params)

# ReferenceFEs
k = params[:k]
T = Float64
model = params[:model]
D = num_cell_dims(model)
reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
reffe_p = ReferenceFE(lagrangian,T,k-1;space=:P)
reffe_j = ReferenceFE(raviart_thomas,T,k-1)
reffe_φ = ReferenceFE(lagrangian,T,k-1)

# Test spaces
Ωf = _fluid_mesh(model,params[:fluid][:domain])
V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
V_p = TestFESpace(Ωf,reffe_p;conformity=p_conformity(Ωf))
V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=BlockMultiFieldStyle())

# Trial spaces

z = zero(VectorValue{D,Float64})
u_bc = params[:bcs][:u][:values]
j_bc = params[:bcs][:j][:values]
U_u = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
U_j = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
U_p = TrialFESpace(V_p)
U_φ = TrialFESpace(V_φ)
U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=BlockMultiFieldStyle())

res, jac = weak_form(params,k)
Tm = params[:matrix_type]
Tv = params[:vector_type]
assem = SparseMatrixAssembler(Tm,Tv,U,V)
op = FEOperator(res,jac,U,V,assem)
xh = zero(U)

r = residual(op,xh)
