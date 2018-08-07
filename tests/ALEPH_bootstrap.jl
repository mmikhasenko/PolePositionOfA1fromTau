# 3
data_bstp = readdlm(joinpath("tau-3pi-analysis","data","aleph_3pi_bstp.txt"));

# data_bstp

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
cdata, ccov = filter_data_and_cov!(data_3pi, cov_3pi, 7, -45)

push!(LOAD_PATH, "src")
using JPACModel
using Fitting

χ2_qtb_disp = build_chi2(intens_qtb_disp, cdata, ccov)
fitres_qtb_disp = fit_with_grad(χ2_qtb_disp, starting_pars=[1.2, 7.45, 0.02], algorithm=:LD_MMA)


vcat(fitres_qtb_disp[1],fitres_qtb_disp[2])

function make_bootstrap(model, data, ccov, bdata, pars=[1.2, 7.45, 0.02])
    _tdata = copy(data)
    ## the bootstrap
    Nb = size(bdata,2)
    results = Vector(Nb);
    @progress for i in 1:Nb
        _tdata[:,2] .= bdata[:,i]
        _cdata, _ccov = filter_data_and_cov!(_tdata, ccov, 7, -45)
        χ2_qtb_disp = build_chi2(model, _cdata, _ccov)
        fitres_qtb_disp = fit_with_grad(χ2_qtb_disp, starting_pars=pars, algorithm=:LD_SLSQP, verbose=false)
        results[i] = vcat(fitres_qtb_disp[1],fitres_qtb_disp[2])
    end
    results
end

# execute the command
@time bresults = make_bootstrap(intens_qtb_disp, tdata, tcov, data_bstp, fitres_qtb_disp[2])

function write_results_to_file(file_name, res)
    bresults_shaped = hcat(res...).'
    writedlm(file_name, bresults_shaped)
    bresults_shaped
end

res_qtb_disp = write_results_to_file("/tmp/qtb_disp_bootstrp.txt", bresults)

using Plots

using SomeRoutines

qtbi = build_function_from_lookup("precalc/qtb_rr5_rect")

function get_qtb_pole(pars, start_xy, qtb_external)
    # find the pole position using interpolated map
    to_minimize(x,y) = abs(1.0/t_qtb_disp_II_external(x+1im*y, pars[1:2]..., qtb_external))
    (v, pars, stat) = fit(to_minimize; starting_pars=start_xy, verbose=false)
    sqrtsp = sqrt(complex(pars...))
end

pole_from_the_main_fit = get_qtb_pole(fitres_qtb_disp[2], [1.3,-0.8], qtbi)
spurious_pole_from_mft = get_qtb_pole(fitres_qtb_disp[2], [0.7,-0.5], qtbi)


poles = Vector{}(size(brest,1))
@progress for i in 1:size(brest,1)
    poles[i] = get_qtb_pole(brest[i,2:3],collect(reim(pole_from_the_main_fit^2)), qtbi)
end

spurs = Vector{}(size(brest,1))
@progress for i in 1:size(brest,1)
    spurs[i] = get_qtb_pole(brest[i,2:3],collect(reim(spurious_pole_from_mft^2)), qtbi)
end

shaped_poles = write_results_to_file("/tmp/qtb_disp_bootstrp_poles.txt", collect.(reim.(poles.^2)))
shaped_spurs = write_results_to_file("/tmp/qtb_disp_bootstrp_spurs.txt", collect.(reim.(spurs.^2)))

readdlm("/tmp/qtb_disp_bootstrp_poles.txt")
readdlm("/tmp/symm_disp_bootstrp_poles.txt")

readdlm("/tmp/symm_disp_bootstrp_poles.txt")

# shaped_spurs
let (cal, grid_x, grid_y) = get_results("qtb_rr5_rect")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
end
scatter!(shaped_poles[:,1], shaped_poles[:,2])
scatter!(shaped_spurs[:,1], shaped_spurs[:,2])
let v1 = pole_from_the_main_fit^2,  v2 = spurious_pole_from_mft^2
    vline!(real.([v1, v2]))
    hline!(imag.([v1, v2]))
end
###################################################################################################
###################################################################################################
###################################################################################################

χ2_symm_disp = build_chi2(intens_symm_disp, cdata, ccov)
fitres_symm_disp = fit_with_grad(χ2_symm_disp, starting_pars=[1.2, 7.45, 0.02], algorithm=:LD_MMA)

@time bsymmresults = make_bootstrap(intens_symm_disp, tdata, tcov, data_bstp, fitres_symm_disp[2])

bsymmresults = [vcat(collect(v)...) for v in bsymmresults]

res_symm_disp = write_results_to_file("/tmp/symm_disp_bootstrp.txt", bsymmresults)

histogram(res_symm_disp[:,1])
# bsymmresults

interfi = build_function_from_lookup("new_interf_rr5")
symmi(s) = qtbi(s)+interfi(s)

function get_symm_pole(pars, start_xy, symm_external)
    # find the pole position using interpolated map
    to_minimize(x,y) = abs(1.0/t_symm_disp_II_external(x+1im*y, pars[1:2]..., symm_external))
    (v, pars, stat) = fit(to_minimize; starting_pars=start_xy, verbose=false)
    sqrtsp = sqrt(complex(pars...))
end


pole_from_the_main_fit_symm = get_symm_pole(fitres_symm_disp[2], [1.5,-0.7], symmi)
spurious_pole_from_mft_symm = get_symm_pole(fitres_symm_disp[2], [0.4,-0.45], symmi)

poles_symm = Vector{}(size(res_symm_disp,1))
@progress for i in 1:size(res_symm_disp,1)
    poles_symm[i] = get_symm_pole(res_symm_disp[i,2:3],collect(reim(pole_from_the_main_fit_symm^2)), symmi)
end

spurs_symm = Vector{}(size(res_symm_disp,1))
@progress for i in 1:size(res_symm_disp,1)
    spurs_symm[i] = get_symm_pole(res_symm_disp[i,2:3],collect(reim(spurious_pole_from_mft_symm^2)), symmi)
end

shaped_poles_symm = write_results_to_file("/tmp/symm_disp_bootstrp_poles.txt", collect.(reim.(poles_symm.^2)))
shaped_spurs_symm = write_results_to_file("/tmp/symm_disp_bootstrp_spurs.txt", collect.(reim.(spurs_symm.^2)))

let (cal, grid_x, grid_y) = get_results("new_interf_rr5")
    contour(grid_x, grid_y, log.(abs.(cal)), levels=50)
end
scatter!(shaped_poles_symm[:,1], shaped_poles_symm[:,2])
scatter!(shaped_spurs_symm[:,1], shaped_spurs_symm[:,2])
let v1 = pole_from_the_main_fit_symm^2,  v2 = spurious_pole_from_mft_symm^2
    vline!(real.([v1, v2]))
    hline!(imag.([v1, v2]))
end
