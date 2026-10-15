include("common.jl")
# Explicit, reproducible rho-radius variation; no editing constants between sessions.
using QuadGK: quadgk
for radius in (SMOKE ? (5.0,) : (3.0,5.0,7.0))
    p0 = sqrt(Masses.mρ2/4-Masses.mπ2)
    function unnormalized(u)
        p = sqrt(u/4-Masses.mπ2)
        ff = p^2/p0^2*(radius^-2+p0^2)/(radius^-2+p^2)
        width = Masses.mρ*Isobars.Γρ*(p/sqrt(u))/(p0/Masses.mρ)*ff
        ff/((Masses.mρ2-u)^2+width^2)*p/(4π*sqrt(u))
    end
    norm = quadgk(u->unnormalized(u)/(2π),4Masses.mπ2,Inf)[1]
    function rho(s)
        s <= 9Masses.mπ2 && return 0.0
        quadgk(u->unnormalized(u)/norm*sqrt(max(0,Isobars.λ(s,u,Masses.mπ2)))/(8π*s)/(2π),
            4Masses.mπ2,(sqrt(s)-Masses.mπ)^2)[1]
    end
    table = SomeRoutines.make_lookup_with_asymptotics_and_tan_map(rho,9Masses.mπ2,1e4)
    model = DispersiveLineshape(table;threshold=9Masses.mπ2,subtraction=(Masses.mπ+Masses.mρ)^2)
    xs = range(0.3,3,length=SMOKE ? 5 : 100)
    write_table("rho_radius_$(radius).csv",["s_GeV2","rho","intensity"],
        permutedims(hcat([[s,table(s),intensity(model,s,1.3,7.0)] for s in xs]...)))
end
