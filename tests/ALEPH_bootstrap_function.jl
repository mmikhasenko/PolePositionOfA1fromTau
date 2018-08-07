using Plots
theme(:juno)

push!(LOAD_PATH, "src")
using JPACModel
using Fitting
using SomeRoutines

using JLD

###################################################################################################
###################################################################################################
###################################################################################################

function make_bootstrap(model, cdata, ccov, bdata, pars=[1.2, 7.45, 0.02])
    _cdata = copy(cdata)
    ## the bootstrap
    Nb = size(bdata,2)
    results = Vector(Nb);
    @progress for i in 1:Nb
        _cdata[:,2] .= bdata[:,i]
        χ2_qtb_disp = build_chi2(model, _cdata, ccov)
        fitres_qtb_disp = fit_with_grad(χ2_qtb_disp, starting_pars=pars, algorithm=:LD_SLSQP, verbose=false)
        results[i] = vcat(collect(fitres_qtb_disp[1:2])...)
    end
    results
end

function write_results_to_file(file_name, res)
    bresults_shaped = hcat(res...).'
    writedlm(file_name, bresults_shaped)
    bresults_shaped
end

function get_pole(model_II, pars, start_xy)
    # find the pole position using interpolated map
    to_minimize(x,y) = abs(1.0/model_II(x+1im*y, pars[1:2]...))
    (v, pars, stat) = fit(to_minimize; starting_pars=start_xy, verbose=false)
    complex(pars...)
end

function bootstrap_chain(model_re, model_II, cdata, ccov, bstpdata, output)
    println("Make the main fit")
    # make the main fit
    χ2 = build_chi2(model_re, cdata, ccov)
    fitres = fit_with_grad(χ2, starting_pars=[1.205, 6.7, 0.017], algorithm=:LD_MMA)

    println("Get Poles from the main fit")
    pole_from_the_main_fit = get_pole(model_II, fitres[2], [1.3,-1.0])
    spurious_pole_from_mft = get_pole(model_II, fitres[2], [0.6,-0.3])
    println("---------> pole_from_the_main_fit = $pole_from_the_main_fit")
    println("---------> spurious_pole_from_mft = $spurious_pole_from_mft")

    println("Make the bootstrap fit")
    @time fits = make_bootstrap(model_re, cdata, ccov, bstpdata, fitres[2])
    fits_shaped = write_results_to_file(output*"_fits.txt", fits)
    # fits_shaped = readdlm(output*"_fits.txt")

    println("Get Poles from the bootstrap")
    println("---------> resonance pole")
    poles = Vector{}(size(fits_shaped,1))
    @progress for i in 1:size(fits_shaped,1)
        poles[i] = get_pole(model_II, fits_shaped[i,2:3], collect(reim(pole_from_the_main_fit)))
    end
    println("---------> spurious pole")
    spurs = Vector{}(size(fits_shaped,1))
    @progress for i in 1:size(fits_shaped,1)
        spurs[i] = get_pole(model_II, fits_shaped[i,2:3], collect(reim(spurious_pole_from_mft)))
    end

    write_results_to_file(output*"_poles.txt", collect.(reim.(poles)))
    write_results_to_file(output*"_spurs.txt", collect.(reim.(spurs)))

    return fitres, pole_from_the_main_fit, spurious_pole_from_mft
end


###################################################################################################
###################################################################################################
###################################################################################################

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

data_bstp = readdlm(joinpath("tau-3pi-analysis","data","aleph_3pi_bstp.txt"));

###################################################################################################

calculate_on_grid(s->qtb(s,5),
            linspace(0.3, 2.0, 100), linspace(-1.2, 0.1, 130),
            "precalc/qtb_rr5_rect")

qtbi = build_function_from_lookup("precalc/qtb_rr5_rect")
model_II_qtb_disp(s, pars...) = t_qtb_disp_II_external(s, pars..., qtbi)

@time rb_qtb_disp = bootstrap_chain(intens_qtb_disp, model_II_qtb_disp, cdata, ccov, data_bstp[7:end-7,:],
    "precalc/qtb_disp_bootstrap")

let sx = linspace(0.3,1.5,20), sy = linspace(-1.2,0.1,20)
    cal = [(y<0 ? model_II_qtb_disp : t_qtb_disp)(x+1im*y,  1.22897, 7.78118) for y in sy, x in sx]
    contour(sx, sy, log.(abs.(cal)))
end
savefig("precalc/qtb_disp.pdf")

##
save("precalc/qtb_disp_bootstrap.jld",
     "main_fit", rb_qtb_disp[1], "pole", rb_qtb_disp[2], "spur", rb_qtb_disp[3])

###################################################################################################

calculate_on_grid(s->interference(s,5),
                  linspace(0.3, 2.0, 100), linspace(-1.3, 0.1, 140),
                  "precalc/new_interf_rr5")

interfi = build_function_from_lookup("precalc/new_interf_rr5")
symmi(s) = qtbi(s)+interfi(s)
model_II_symm_disp(s, pars...) = t_symm_disp_II_external(s, pars..., symmi)

@time rb_symm_disp = bootstrap_chain(intens_symm_disp, model_II_symm_disp, cdata, ccov, data_bstp[7:end-7,:],
    "precalc/symm_disp_bootstrap")

@time let sx = linspace(0.3,1.5,20), sy = linspace(-0.9,0.1,20)
    cal = [(y<0 ? model_II_symm_disp : t_symm_disp)(x+1im*y,  1.1949991980985653, 5.990846297030663) for y in sy, x in sx]
    contour(sx, sy, log.(abs.(cal)))
end
savefig("precalc/symm_disp.pdf")

###################################################################################################

save("precalc/symm_disp_bootstrap.jld",
     "main_fit", rb_symm_disp[1], "pole", rb_symm_disp[2], "spur", rb_symm_disp[3])

###################################################################################################
