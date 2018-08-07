# complete data
data_3pi = readdlm(joinpath("tau-3pi-analysis","data","ALEPH.2005.3pi.formatted"));
corr_3pi = readdlm(joinpath("tau-3pi-analysis","data","stat_corr_3pi.txt"));
cov_3pi = let n=size(data_3pi,1)
    [data_3pi[i,3]*corr_3pi[i,j]/100*data_3pi[j,3] for i in 1:n, j in 1:n]
end

function filter_data_and_cov!(data, mat, first_accepted=1, last_accepted=0)
    # first map: nonzero
    trfilt = [v!=0.0 for v in data[:,2]]
    # second map: truncation
    cfilt = first_accepted:(sum(trfilt)+last_accepted)

    # application
    return data[trfilt,:][cfilt,:], mat[trfilt,trfilt][cfilt,cfilt]
end
tdata, tcov = filter_data_and_cov!(data_3pi, cov_3pi)
cdata, ccov = filter_data_and_cov!(data_3pi, cov_3pi, 7, -7)

push!(LOAD_PATH, "src")
using JPACModel
using Fitting

χ2_qtb_disp = build_chi2(intens_qtb_disp, cdata, ccov)
fitres_qtb_disp = fit_with_grad(χ2_qtb_disp, starting_pars=[1.2, 7.0, 0.02])

χ2_symm_disp = build_chi2(intens_symm_disp, cdata, ccov)
fitres_symm_disp = fit_with_grad(χ2_symm_disp, starting_pars=[1.2, 7.0, 0.02])

using JLD
save("tau-3pi-analysis/syst_test_rho_shape/G-40.0.jld", "qtb_disp", fitres_qtb_disp, "symm_disp", fitres_symm_disp)


# ######################################################################
#
# Rv = [1,2,2.5,3,3.5,4,5,6,7,8,10]
# χQTBDISP = Vector{Float64}(0)
# χSYMMDISP = Vector{Float64}(0)
# for R in Rv
#     d = load("tau-3pi-analysis/syst_test_rho_shape/R$(R).jld");
#     fitres_qtb_disp = d["qtb_disp"]
#     fitres_symm_disp = d["symm_disp"]
#     # @show fitres_qtb_disp[1]
#     push!(χQTBDISP, fitres_qtb_disp[1])
#     push!(χSYMMDISP, fitres_symm_disp[1])
# end
#
# using Plots
#
# plot(Rv, χQTBDISP, lab="QTB-DISP", xlab="R (1/GeV)", ylab="chi2")
# plot!(Rv, χSYMMDISP, lab="SYMM-DISP")
# savefig("tau-3pi-analysis/syst_test_rho_shape/R-summary.png")
#
# ######################################################################
#
# Gv = [-40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0];
# χQTBDISP = Vector{Float64}(0)
# χSYMMDISP = Vector{Float64}(0)
# for G in Gv
#     d = load("tau-3pi-analysis/syst_test_rho_shape/"*(
#             G ≥ 0.0 ? "G+$(abs(G))" : "G-$(abs(G))"
#         )*".jld");
#     fitres_qtb_disp = d["qtb_disp"]
#     fitres_symm_disp = d["symm_disp"]
#     push!(χQTBDISP, fitres_qtb_disp[1])
#     push!(χSYMMDISP, fitres_symm_disp[1])
# end
#
# using Plots
#
# plot(Gv, χQTBDISP, lab="QTB-DISP", xlab="shift to the rho width (MeV)", ylab="chi2")
# plot!(Gv, χSYMMDISP, lab="SYMM-DISP")
# savefig("tau-3pi-analysis/syst_test_rho_shape/G-summary.png")
#
# #######################################################################
#
# Mv = [-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0];
# χQTBDISP = Vector{Float64}(0)
# χSYMMDISP = Vector{Float64}(0)
# for M in Mv
#     d = load("tau-3pi-analysis/syst_test_rho_shape/"*(
#             M ≥ 0.0 ? "M+$(abs(M))" : "M-$(abs(M))"
#         )*".jld");
#     fitres_qtb_disp = d["qtb_disp"]
#     fitres_symm_disp = d["symm_disp"]
#     push!(χQTBDISP, fitres_qtb_disp[1])
#     push!(χSYMMDISP, fitres_symm_disp[1])
# end
#
# using Plots
#
# plot(Mv, χQTBDISP, lab="QTB-DISP", xlab="shift to the rho mass (MeV)", ylab="chi2")
# plot!(Mv, χSYMMDISP, lab="SYMM-DISP")
# savefig("tau-3pi-analysis/syst_test_rho_shape/M-summary.png")
