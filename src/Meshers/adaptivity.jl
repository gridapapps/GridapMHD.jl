
@enum AdaptivityMethod begin
  NONE = 0
  SIMPLEXIFY = 1
  BARYCENTER = 2
  SIMPLEXIFY_AND_BARYCENTER = 3
end

function adapt_mesh(model,method)
  method = AdaptivityMethod(method)
  if method == SIMPLEXIFY
    model = simplexify(model)
  elseif method == BARYCENTER
    model = Gridap.Adaptivity.get_model(
      Gridap.Adaptivity.refine(model, refinement_method = "barycentric")
    )
  elseif method == SIMPLEXIFY_AND_BARYCENTER
    model = simplexify(model)
    model = Gridap.Adaptivity.get_model(
      Gridap.Adaptivity.refine(model, refinement_method = "barycentric")
    )
  end
  return model
end
