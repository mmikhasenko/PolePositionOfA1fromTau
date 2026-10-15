using PolePositionOfA1fromTau
using PolePositionOfA1fromTau: JPACModel, Masses, Isobars, Fitting, SomeRoutines
using DelimitedFiles, LinearAlgebra, Random, TOML
const ROOT = dirname(@__DIR__)
const OUT = get(ENV, "A1_OUTPUT_DIR", joinpath(ROOT, "artifacts"))
const SMOKE = "--smoke" in ARGS
mkpath(OUT)

function write_table(name, header, rows)
    all(isfinite, rows) || error("Non-finite values in output $name")
    open(joinpath(OUT,name), "w") do io
        println(io, join(header, ','))
        writedlm(io, rows, ',')
    end
end

function aleph_data()
    data = readdlm(joinpath(ROOT,"data","ALEPH.2005.3pi.formatted"))
    corr = readdlm(joinpath(ROOT,"data","stat_corr_3pi.txt"))
    cov = Diagonal(data[:,3])*corr*Diagonal(data[:,3])/100
    nonzero = findall(!iszero, data[:,2])
    selected = nonzero[7:end-7]
    data[selected,:], cov[selected,selected], selected
end

# Precompute only quantities independent of fit parameters; preserve AD in m,g,n.
function prepared_fit(kind=:qtb; dispersive=true)
    data, cov, indices = aleph_data()
    rho = JPACModel.density(kind)
    model = JPACModel.dispersive_model(kind)
    xs = data[:,1]
    loops = dispersive ? selfenergy.(Ref(model),xs) : im .* rho.(xs)
    phase = rho.(xs) .* (Masses.mτ^2 .- xs).^2 ./ xs .* (1 .+ 2xs/Masses.mτ^2)
    precision = inv(Symmetric(cov))
    function objective(y)
        function chi2(m,g,n)
            predicted = abs2.(n*g^2 ./ (m^2 .- xs .- g^2 .* loops ./ 2)) .* phase
            residual = y-predicted
            dot(residual,precision*residual)
        end
    end
    objective, data, indices
end
function fit_aleph(kind=:qtb; dispersive=true)
    objective, data, indices = prepared_fit(kind; dispersive)
    initial = Fitting.fit(objective(data[:,2]); starting_pars=[1.2,7.0,0.02], algorithm=:LN_NELDERMEAD, verbose=false)
    result = Fitting.fit_with_grad(objective(data[:,2]); starting_pars=initial[2], algorithm=:LD_LBFGS, verbose=false)
    label = string(kind, dispersive ? "_disp" : "")
    open(joinpath(OUT,"fit_"*label*".toml"),"w") do io
        TOML.print(io, Dict("chi2"=>result[1],"parameters"=>result[2],
            "parameter_names"=>["m_GeV","g_GeV","normalization"],"status"=>string(result[3]),
            "fit_rows"=>indices,"covariance"=>[collect(row) for row in eachrow(result[4])]))
    end
    println(label, ": chi2=",result[1]," parameters=",result[2]," status=",result[3])
    result
end
function export_lineshape(; kind=:qtb, n=SMOKE ? 8 : 150)
    model = JPACModel.dispersive_model(kind)
    xs = range(9Masses.mπ2+0.001,3.0,length=n)
    rows = [begin
        loop=selfenergy(model,s); t=amplitude(model,s,1.3,7.0)
        [s,model.rho(s),real(loop),imag(loop),real(t),imag(t),intensity(model,s,1.3,7.0)]
    end for s in xs]
    write_table("lineshape_$(kind).csv",["s_GeV2","rho","Pi_re","Pi_im","t_re","t_im","intensity"],permutedims(hcat(rows...)))
    open(joinpath(OUT,"lineshape_$(kind).toml"),"w") do io
        TOML.print(io,Dict("model"=>string(kind),"m_GeV"=>1.3,"g_GeV"=>7.0,
            "normalization"=>1.0,"parameters_are"=>"illustrative, not a refit",
            "threshold_GeV2"=>model.threshold,"subtraction_GeV2"=>model.subtraction,
            "density_grid_points"=>100,"density_high_tail"=>1/(8π),
            "density_tail_starts_GeV2"=>1e4,"julia_version"=>string(VERSION)))
    end
end
function export_grid(; kind=:qtb)
    f = kind == :qtb ? qtb : interference
    n = SMOKE ? 3 : 50
    SomeRoutines.calculate_on_grid(f,range(0.3,2.0,length=n),range(-0.8,-0.01,length=n),joinpath(OUT,"$(kind)_complex"))
end
