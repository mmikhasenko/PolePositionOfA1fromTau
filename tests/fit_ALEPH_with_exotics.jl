# complete data
data_3pi = readdlm(joinpath("data","ALEPH.2005.3pi.formatted"));
corr_3pi = readdlm(joinpath("data","stat_corr_3pi.txt"));
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

using sQTBDISPk

χ2_sqtb_dispk = build_chi2(intens_sqtb_dispk, cdata, ccov)
fitres_sqtb_dispk2 = fit_with_grad(χ2_sqtb_dispk;
    starting_pars=[1.8, 16.0, 0.017],
    algorithm=:LD_MMA)

fitres_sqtb_dispk3 = fit_with_grad(χ2_sqtb_dispk, starting_pars=[1.23, 7.4, -0.03, 0.0162],
    algorithm=:LD_MMA)

using Plots
theme(:juno)

function plot_vs_data(pars)
    plot(cdata[:,1], cdata[:,2], yerr=cdata[:,3], l=nothing)
    plot!(x->intens_sqtb_dispk(x,pars...), linspace(cdata[1,1],cdata[end,1],100))
end

function plot_vs_the_second_sheet(pars)
    sxv = linspace(0.3,2.0,30)
    syv = linspace(-1.0, 0.1, 30)
    cal = [(sy>=0 ? t_sqtb_dispk(sx+1im*sy, pars[1:(end-1)]...) :
                   t_sqtb_dispk_II(sx+1im*sy, pars[1:(end-1)]...)) for sy in syv, sx in sxv]
    contour(sxv, syv, log.(abs.(cal)), levels=100, colorbar=false)
end

plot_vs_data(fitres_sqtb_dispk3[2])
plot_vs_the_second_sheet(fitres_sqtb_dispk3[2])

fitres_sqtb_dispk4 = fit_with_grad(χ2_sqtb_dispk, starting_pars=[1.23, 6.0, 0.0, -40, 0.016],
    algorithm=:LD_SLSQP)

plot_vs_data(fitres_sqtb_dispk4[2])
plot_vs_the_second_sheet(fitres_sqtb_dispk4[2])

using SomeRoutines

qtb_lookup = build_function_from_lookup("qtb_rr5_rect")
function plot_vs_the_second_sheet_external(pars)
    sxv = linspace( 0.3, 2.0, 100)
    syv = linspace(-1.0, 0.1, 130)
    cal = [(sy>=0 ? t_sqtb_dispk(sx+1im*sy, pars[1:(end-1)]...) :
                    t_sqtb_dispk_II_external(qtb_lookup, sx+1im*sy, pars[1:(end-1)]...)) for sy in syv, sx in sxv]
    contour(sxv, syv, log.(abs.(cal)), levels=100, colorbar=false)
end

using JLD
function save_the_second_sheet_external(pars,output)
    sxv = linspace( 0.3, 2.0, 100)
    syv = linspace(-1.2, 0.1, 130)
    cal = [(sy>=0 ? t_sqtb_dispk(sx+1im*sy, pars[1:(end-1)]...) :
                    t_sqtb_dispk_II_external(qtb_lookup, sx+1im*sy, pars[1:(end-1)]...)) for sy in syv, sx in sxv]
    save(output, "sx", sxv, "sy", syv, "amp", cal)
end


@time save_the_second_sheet_external(fitres_sqtb_dispk2[2], "/tmp/sqtb_dispk2.jld")
@time save_the_second_sheet_external(fitres_sqtb_dispk3[2], "/tmp/sqtb_dispk3.jld")
@time save_the_second_sheet_external(fitres_sqtb_dispk4[2], "/tmp/sqtb_dispk4.jld")
