
function gmg_solver(mh,trials,tests,biform,qdegree)
  smatrices = compute_gmg_matrices(mh,trials,tests,biform,qdegree)
  restrictions, prolongations = setup_transfer_operators(trials,
                                                         qdegree;
                                                         mode=:residual,
                                                         solver=IS_ConjugateGradientSolver(;reltol=1.e-6))

  smoothers = gmg_patch_smoothers(mh,tests,biform,qdegree)

  # TODO: Reuse top level matrix in gmg. In fact, all matrices should inputed (and modified) by numerical_setup
  gmg = GMGLinearSolver(mh,
                        smatrices,
                        prolongations,
                        restrictions,
                        pre_smoothers=smoothers,
                        post_smoothers=smoothers,
                        coarsest_solver=PETScLinearSolver(petsc_mumps_setup),
                        maxiter=1,
                        rtol=1.0e-8,
                        verbose=false,
                        mode=:preconditioner)
  solver = GMRESSolver(5;Pr=gmg,m_add=3,maxiter=15,rtol=1.0e-8,verbose=i_am_main(parts))
  return solver
end

function compute_gmg_matrices(mh,trials,tests,biform,qdegree)
  nlevs = num_levels(trials)

  mats = Vector{PSparseMatrix}(undef,nlevs)
  for lev in 2:nlevs
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      model = get_model(mh,lev)
      U = get_fe_space(trials,lev)
      V = get_fe_space(tests,lev)
      Ω = Triangulation(model)
      dΩ = Measure(Ω,qdegree[lev])
      a(u,v) = biform(u,v,dΩ)
      mats[lev] = assemble_matrix(ai,U,V)
    end
  end
  return mats
end

function gmg_patch_smoothers(mh,tests,biform,qdegree)
  patch_decompositions = PatchDecomposition(mh)
  patch_spaces = PatchFESpace(tests,patch_decompositions)

  nlevs = num_levels(mh)
  smoothers = Vector{RichardsonSmoother}(undef,nlevs-1)
  for lev in 1:nlevs-1
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      PD = patch_decompositions[lev]
      Ph = get_fe_space(patch_spaces,lev)
      Vh = get_fe_space(tests,lev)
      Ω  = Triangulation(PD)
      dΩ = Measure(Ω,qdegree[lev])
      local_solver   = LUSolver()
      patch_smoother = PatchBasedLinearSolver(biform,Ph,Vh,dΩ,local_solver)
      smoothers[lev] = RichardsonSmoother(patch_smoother,5,0.2)
    end
  end
  return smoothers
end
