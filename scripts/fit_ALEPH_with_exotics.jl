include("common.jl")
using PolePositionOfA1fromTau: sQTBDISPk
# Fit the three K-matrix variants defined in src/sQTBDISPk.jl.
data,cov,_ = aleph_data()
for start in (SMOKE ? ([1.2,7.0,0.02],) : ([1.2,7.0,0.02],[1.2,7.0,0.1,0.02],[1.2,7.0,2.5,0.1,0.02]))
    chi2 = Fitting.build_chi2(sQTBDISPk.intens_sqtb_dispk,data,cov)
    result = Fitting.fit(chi2; starting_pars=start,algorithm=:LN_NELDERMEAD,verbose=false)
    open(joinpath(OUT,"fit_kmatrix_$(length(start)-1).toml"),"w") do io
        TOML.print(io,Dict("chi2"=>result[1],"parameters"=>result[2],"status"=>string(result[3])))
    end
end
