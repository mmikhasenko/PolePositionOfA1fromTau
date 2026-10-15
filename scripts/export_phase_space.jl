include("common.jl")
xs = range(9Masses.mπ2,3.0,length=SMOKE ? 8 : 150)
rows = [[s,qtb(s),symm(s),s > (Masses.mρ+Masses.mπ)^2 ? sqrt(Isobars.λ(s,Masses.mρ2,Masses.mπ2))/(8π*s) : 0.0] for s in xs]
write_table("phase_space.csv",["s_GeV2","QTB","SYMM","stable_rho_pi"],permutedims(hcat(rows...)))
