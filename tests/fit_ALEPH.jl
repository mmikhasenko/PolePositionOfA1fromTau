# complete data
# data_3pi = readdlm("C:\\Users\\mikha\\cernbox\\current_projects\\unitarity\\isospin\\precalc\\ALEPH.2005.3pi.formatted");
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


χ2_qtb = build_chi2(intens_qtb, cdata, ccov)
fitres_qtb = fit_with_grad(χ2_qtb, starting_pars=[1.2, 7.0, 0.6])

χ2_qtb_disp = build_chi2(intens_qtb_disp, cdata, ccov)
fitres_qtb_disp = fit_with_grad(χ2_qtb_disp, starting_pars=[1.2, 7.0, 0.6])

χ2_symm = build_chi2(intens_symm, cdata, ccov)
fitres_symm = fit_with_grad(χ2_symm, starting_pars=[1.2, 7.0, 0.6])

χ2_symm_disp = build_chi2(intens_symm_disp, cdata, ccov)
fitres_symm_disp = fit_with_grad(χ2_symm_disp, starting_pars=[1.2, 7.0, 0.6])

fitres_qtb_disp

using Plots
# theme()
gr()
plot(tdata[:,1], tdata[:,2], yerr=tdata[:,3], l=nothing, lab="")
plot!(cdata[:,1], cdata[:,2], yerr=cdata[:,3], l=nothing, lab="ALEPH")

let sv = linspace(tdata[1,1],tdata[end,1],100), pars = fitres_symm[2]
    cal = intens_symm.(sv, pars...)
    range = [cdata[1,1]<v<tdata[end,1] for v in sv]
    plot!(sv, cal, lab="", l=(:dash)) # full range
    plot!(sv[range], cal[range], lab="symm") # fit range
end

let sv = linspace(tdata[1,1],tdata[end,1],100), pars = fitres_symm_disp[2]
    cal = intens_symm_disp.(sv, pars...)
    range = [cdata[1,1]<v<tdata[end,1] for v in sv]
    plot!(sv, cal, lab="", l=(:dash)) # full range
    plot!(sv[range], cal[range], lab="symm_disp") # fit range
end

let sv = linspace(tdata[1,1],tdata[end,1],100), pars = fitres_qtb_disp[2]
    cal = intens_qtb_disp.(sv, pars...)
    range = [cdata[1,1]<v<tdata[end,1] for v in sv]
    plot!(sv, cal, lab="", l=(:dash)) # full range
    plot!(sv[range], cal[range], lab="qtb_disp") # fit range
end


# plot the second sheet
@time let sx = 1.3, syv = linspace(-1.0, 0.1, 99), pars = fitres_symm[2]
    cal = [t_symm(sx+1im*sy, pars[1:2]...) for sy in syv]
    plot!(syv, real(cal), lab="") # full range
    plot!(syv, imag(cal), lab="") # full range
end

@time let sx = 1.3, syv = linspace(-1.0, 0.1, 99), pars = fitres_symm_disp[2]
    cal = [(sy>0 ? t_symm_disp : t_symm_disp_II)(sx+1im*sy, pars[1:2]...) for sy in syv]
    plot(syv, real(cal), lab="") # full range
    plot!(syv, imag(cal), lab="") # full range
end

@time let sxv = linspace(0.3,2.0,200), syv = linspace(-1.0, 0.1, 200), pars = fitres_symm_disp[2]
    cal = [(sy>0 ? t_symm_disp(sx+1im*sy, pars[1:2]...):
                   t_symm_disp_II_external(sx+1im*sy, pars[1:2]...,
                        s->(qtbi(s)+interfi(s)))) for sy in syv, sx in sxv]
    contour(sxv, syv, log.(abs.(cal)), levels=100, colorbar=false)
end

t_symm_disp_II_external(1.1-0.6im, fitres_symm_disp[2][1:2]...,s->(qtbi(s)+interfi(s)))
t_symm_disp_II(1.1-0.6im, fitres_symm_disp[2][1:2]...)

using Plots
theme(:juno)

@time let sxv = linspace(0.3,2.0,10), syv = linspace(-1.0, 0.1, 100), pars = fitres_qtb_disp[2]
    cal = [(sy>=0 ? t_qtb_disp(sx+1im*sy, pars[1:2]...) :
                   t_qtb_disp_II(sx+1im*sy, pars[1:2]...)) for sy in syv, sx in sxv]
    contour(sxv, syv, log.(abs.(cal)), levels=100, colorbar=false)
end

# find the pole position using interpolated map
to_minimize(x,y) = abs(1.0/t_symm_disp_II_external(x+1im*y, fitres_symm_disp[2][1:2]...,s->(qtbi(s)+interfi(s))))
let (v, pars, stat) = fit(to_minimize; starting_pars=[1.5,-0.75])
    sqrtsp = sqrt(complex(pars...))
    println("pole mass is $(real(sqrtsp)), pole width is $(-2*imag(sqrtsp))")
end # pole mass is 1.210701837808158, pole width is 0.582255621800543

# find the pole position using exact function
to_minimize(x,y) = abs(1.0/t_symm_disp_II(x+1im*y, fitres_symm_disp[2][1:2]...))
let (v, pars, stat) = fit(to_minimize; starting_pars=[1.5,-0.75])
    sqrtsp = sqrt(complex(pars...))
    println("pole mass is $(real(sqrtsp)), pole width is $(-2*imag(sqrtsp))")
end # pole mass is 1.2106825505711496, pole width is 0.5822446532613349
# --------------------------------------------------------------------------



# find the pole position using exact function
to_minimize(x,y) = abs(1.0/t_qtb_disp_II(x+1im*y, fitres_qtb_disp[2][1:2]...))
let (v, pars, stat) = fit(to_minimize; starting_pars=[1.5,-0.75])
    sqrtsp = sqrt(complex(pars...))
    println("pole mass is $(real(sqrtsp)), pole width is $(-2*imag(sqrtsp))")
end # pole mass is 1.166468907756505, pole width is 0.798660853469603
# --------------------------------------------------------------------------

# plot the second sheet
@time let sx = 1.1, syv = linspace(-1.2, 0.1, 99), pars = fitres_qtb_disp[2]
    cal = [(sy>0 ? t_qtb_disp : t_qtb_disp_II)(sx+1im*sy, pars[1:2]...) for sy in syv]
    plot(syv, real(cavbl), lab="") # full range
    plot!(syv, imag(cal), lab="") # full range
end

@time let sxv = linspace(0.7,1.5,50), syv = linspace(-1.1, 0.1, 50), pars = fitres_qtb_disp[2]
    cal = [(sy>0 ? t_qtb_disp : t_qtb_disp_II)(sx+1im*sy, pars[1:2]...) for sy in syv, sx in sxv]
    contour(sxv, syv, log.(abs.(cal)), levels=100, colorbar=false)
end
