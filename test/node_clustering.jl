  using LinearAlgebra
  using Plots
  using Roots

  nx = 40
  ny = 40
  Lx = 1.0
  Ly = 1.0

  O = (0.0,0.0)
  L = (Lx,Ly)
  n = (nx,ny)


  # Choose a nice delta
  # E.g., 1/4 of points in 1/10 of length
  _δ = 1/10 * L[1]
  δ = (_δ, _δ)
  nδ = (10,10) #Int.(floor.(n.*0.25))

  function compute_map(O,L,n,δ,nδ)

    # Reference space [0,2]
    Lr = 2
    h = inv.(n).*Lr

    g(x,α) = ((exp(x)-1)/(exp(1)-1))^α
    r(i) = α -> g(nδ[i]*h[i],α) - δ[i]*Lr
    α = [ find_zero(r(i),1) for i in 1:length(L) ]

    function _geo_map(xi,αi,Li)
      xni = xi
      _xi = xi*2/Li
      _xni = _xi
      f(xi) = αi < 1 ? xi : g(xi,αi)
      xni = _xi <= 1.0 ? f(_xi) : 2 - f(2-_xi)
      xni*Li/2.0
    end

    geo_map = x-> map(_geo_map,x,α,L)
  end

  geo_map = compute_map(O,L,n,δ,nδ)

  hp = L./n

  x = Vector{Float64}[]
  for j in 0:n[2]
    x = vcat(x,[[hp[1]*i,hp[2]*j] for i in 0:n[1]])
  end

  nx = geo_map.(x)

  function plot_points(points::Vector{Vector{Float64}})
    n_dims = length(points[1])
    if n_dims == 2
      plot([p[1] for p in points], [p[2] for p in points], seriestype=:scatter, legend=false)
    elseif n_dims == 3
      plot3d([p[1] for p in points], [p[2] for p in points], [p[3] for p in points], seriestype=:scatter, legend=false)
    else
      error("Plotting is only supported for 2D and 3D points")
    end
  end

  plot_points(x)
  plot_points(nx)
