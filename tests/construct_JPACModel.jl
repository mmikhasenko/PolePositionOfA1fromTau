push!(LOAD_PATH, "src")
using JPACModel
# using Masses

using Plots
let sv = linspace(0.25,3.0,100)
    cal = [8π*qtb(s) for s in sv]
    plot(sv, cal, lab="qtb", xlab="s (GeV^2)")
    cal = [8π*symm(s) for s in sv]
    plot!(sv, cal, lab="symm")
end

let ev = linspace(0.5,2.5,100)
    cal = [8π*qtb(e^2) for e in ev]
    plot(ev, cal, lab="qtb", xlab="e (GeV)")
    cal = [8π*symm(e^2) for e in ev]
    plot!(ev, cal, lab="symm")
end

qtb(1.1+0.1im)
qtb(1.1-0.1im)

Juno.Profile.clear()
@profile interference(1.1+0.1im)
Juno.profiler()

interference(1.1-0.1im)
