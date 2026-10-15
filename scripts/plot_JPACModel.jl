include("common.jl")
ENV["GKSwstype"] = "100"
using Plots
export_lineshape()
data = readdlm(joinpath(OUT,"lineshape_qtb.csv"),',',Float64; skipstart=1)
p = plot(data[:,1],data[:,end]; xlabel="s [GeV²]",ylabel="Intensity [arbitrary units]",label="QTB dispersive (illustrative parameters)")
savefig(p,joinpath(OUT,"lineshape_qtb.pdf"))
