### test
push!(LOAD_PATH, "src")
using JPACModel
using Plots
theme(:juno)
using Masses
using DalitzPlotAnalysis

phase_space_factors = []
sval = linspace(9mπ2, 3.0, 300)
@time let sv = sval
    cal = [s < 9mπ2 ? 0.0 : qtb(s) for s in sv];
    push!(phase_space_factors, cal);
    plot(sv, cal)
    #
    cal = [s < 9mπ2 ? 0.0 : symm(s) for s in sv];
    plot!(sv, cal)
    push!(phase_space_factors, cal);
    #
    cal = [s < (mρ+mπ)^2 ? 0.0 : sqrt(λ(s,mρ^2,mπ2))/(8π*s) for s in sv]
    plot!(sv, cal)
    push!(phase_space_factors, cal);
end

using JLD
save("tau-3pi-analysis/final_plots/phase_space_factors.jld",
    "sv", sval,
    "QTB", phase_space_factors[1],
    "SYMM", phase_space_factors[2],
    "stable", phase_space_factors[3]);
