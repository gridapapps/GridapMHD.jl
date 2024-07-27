using DrWatson
using DataFrames
using Plots
using LinearAlgebra

function log_slope(x,y)
  x,y = log10.(x), log10.(y)
  x = [x ones(length(x))]
  A = (inv(x'*x)*x') * y
  A[1]
end

function plot_error!(x,y;
  name="u",norm="L^2",slope=true,slope_range=1:length(x))

  slope_str = ""
  if slope
    a = log_slope(x[slope_range],y[slope_range])
    slope_str = " \\quad  [$(round(a,digits=2))]"
  end
  label = "\$|| $name - $(name)_h || _{$norm}  $slope_str\$"
  plot!(x,y,label=label,marker=:auto,linestyle=:dash,linewidth=1.5)
end

raw = collect_results(datadir())

r = copy(raw)
r = filter(i-> startswith(i[:title],"transient_") ,r)
r = filter(i->contains(i[:title],"nonlineartime_fespace"),r)
sort!(r,:Δt)

Δt = r[:,:Δt]
euh = map(last,r[:,:uh_el2])
ejh = map(last,r[:,:jh_el2])
eph = map(last,r[:,:ph_el2])
eφh = map(last,r[:,:φh_el2])
euh_h1 = map(last,r[:,:uh_eh1])


slope_range = 1:4
plot()
plot_error!(Δt,euh;name="u",slope_range)
plot_error!(Δt,euh_h1;name="u",norm="H^1",slope_range)
plot_error!(Δt,ejh;name="j",slope_range)
# plot_error!(Δt,eph;name="p",slope_range)
plot_error!(Δt,eφh;name="φ",slope_range)

tiks = exp10.(-15:0)
tiks5 = exp10.(-15:0.5:0)
tiks2 = exp2.(-10:1:0), map(i->"\$2^{$i}\$",(-10:1:0))
tiks2 = exp2.(-10:1:0), map(i->"1/$(2^-i)",(-10:1:0))
plot!(xscale=:log10,yscale=:log10)
plot!(xlabel="Time step size \$(\\Delta t)\$",ylabel="Error Norm")
plot!(xticks=tiks2,yticks=tiks)
plot!(legend=:bottomright)

# plot!(Δt,Δt*1e-2,label="\$\\Delta t ^1 \$",color=:black,linestyle=:dash)
plot!(Δt,Δt.^2*1e-2,label="\$\\Delta t ^2\$",color=:black,linestyle=:dashdot)

savefig("t_convergence.pdf")


# Spatial Convergence
r = copy(raw)
r = filter(i->startswith(i[:title],"transient_") ,r)
r = filter(i->contains(i[:title],"stationary_nonfespace"),r)

sort!(r,:ncells)

n = r[:,:ncells]
h = n .^ (-1/3)


f = last
euh = map(f,r[:,:uh_el2])
ejh = map(f,r[:,:jh_el2])
eph = map(f,r[:,:ph_el2])
eφh = map(f,r[:,:φh_el2])
euh_h1 = map(f,r[:,:uh_eh1])

slope_range = length(h)-1:length(h)
plot()
plot_error!(h,euh;name="u",slope_range)
plot_error!(h,euh_h1;name="u",norm="H^1",slope_range)
plot_error!(h,ejh;name="j",slope_range)
# plot_error!(h,eph;name="p",slope_range)
plot_error!(h,eφh;name="φ",slope_range)

tiks10 = exp10.(-15:10)
tiks2 = exp2.(-10:1:0), map(i->"1/$(2^-i)",(-10:1:0))
plot!(xscale=:log10,yscale=:log10)
plot!(xlabel="Cell size (h)",ylabel="Error Norm")
plot!(xticks=tiks2,yticks=tiks10)
plot!(legend=:bottomright)
plot!(h,h.^2*1e-1,label="\$h^2\$",color=:black,linestyle=:dash)
plot!(h,h.^3*1e-2,label="\$h^3\$",color=:black,linestyle=:dashdot)

savefig("h_convergence_stationary.pdf")


# Spatial Convergence
r = copy(raw)
r = filter(i->startswith(i[:title],"transient_") ,r)
r = filter(i->contains(i[:title],"lineartime_nonfespace"),r)

sort!(r,:ncells)
n = r[:,:ncells]
h = n .^ (-1/3)



f = first
euh = map(f,r[:,:uh_el2])
ejh = map(f,r[:,:jh_el2])
eph = map(f,r[:,:ph_el2])
eφh = map(f,r[:,:φh_el2])
euh_h1 = map(f,r[:,:uh_eh1])

slope_range = 2:length(h)
plot()
plot_error!(h,euh;name="u",slope_range)
plot_error!(h,euh_h1;name="u",norm="H^1",slope_range)
plot_error!(h,ejh;name="j",slope_range)
# plot_error!(h,eph;name="p",slope_range)
plot_error!(h,eφh;name="φ",slope_range)

tiks10 = exp10.(-15:10)
tiks2 = exp2.(-10:1:0), map(i->"1/$(2^-i)",(-10:1:0))
plot!(xscale=:log10,yscale=:log10)
plot!(xlabel="Cell size (h)",ylabel="Error Norm")
plot!(xticks=tiks2,yticks=tiks10)
plot!(legend=:bottomright)
plot!(h,h.^2*1e-1,label="\$h^2\$",color=:black,linestyle=:dash)
plot!(h,h.^3*1e-2,label="\$h^3\$",color=:black,linestyle=:dashdot)

savefig("h_convergence_transient.pdf")


#
# Test: p = sin(x[1]) + offset (zeromean)
# Convergencia temps (NS en gridap)
# Write report, solutions, plots wo
# test Measure k*2*D
