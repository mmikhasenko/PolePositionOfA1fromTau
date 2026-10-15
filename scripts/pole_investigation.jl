include("common.jl")
# This diagnoses the historical QTB continuation, not a certified new pole result.
model = JPACModel.dispersive_model(:qtb)
objective(x,y) = abs2(inv(amplitude(model,complex(x,y),1.3,7.0;
    sheet=:II,continued_density=qtb)))
result = Fitting.fit(objective; starting_pars=[1.3,-0.6],verbose=false)
open(joinpath(OUT,"pole_diagnostic.toml"),"w") do io
    TOML.print(io,Dict("s_re"=>result[2][1],"s_im"=>result[2][2],
        "inverse_amplitude_squared"=>result[1],"status"=>string(result[3]),
        "scope"=>"historical QTB contour; illustrative parameters; global sheet transport not certified"))
end
println(result)
