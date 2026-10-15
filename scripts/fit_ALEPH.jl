include("common.jl")
for kind in (SMOKE ? (:qtb,) : (:qtb,:symm)), dispersive in (SMOKE ? (true,) : (false,true))
    fit_aleph(kind; dispersive)
end
