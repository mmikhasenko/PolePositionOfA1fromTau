# Separate Julia processes catch hidden interactive-session dependencies.
const root=dirname(@__DIR__)
for script in ("fit_example", "export_lineshape", "construct_JPACModel", "export_phase_space",
               "fit_ALEPH", "ALEPH_bootstrap", "precalculate_qtb_interference",
               "isobar_analytical_structure", "plot_JPACModel", "pole_investigation",
               "fit_ALEPH_with_exotics", "systematic_studies_of_rho_shape_influence")
    println("Smoke: ", script); flush(stdout)
    run(`$(Base.julia_cmd()) --project=$root $(joinpath(root,"scripts",script*".jl")) --smoke`)
end
