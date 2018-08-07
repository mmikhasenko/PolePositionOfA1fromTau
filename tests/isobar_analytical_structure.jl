push!(LOAD_PATH, "src")
using Isobars
using Plots
theme(:juno)

gr()
let xv = linspace(-0.5,1,300), yv = linspace(-0.5,0.5,300)
    plot(layout=grid(2,1))
    cal = [f1_I(x+1im*y) for y in yv, x in xv]
    plot!(xv, yv, log.(abs.(cal)), levels=40, subplot=1, colorbar=false)
    cal = [f1_II(x+1im*y) for y in yv, x in xv]
    plot!(xv, yv, log.(abs.(cal)), levels=40, subplot=2, colorbar=false)
end


let xv = linspace(-0.5,1,300), y = 1e-5im
    cal = [f1_I(x+1im*y) for x in xv]
    plot( xv, real.(cal))
    plot!(xv, imag.(cal))
    plot!(xv,  abs.(cal))
end
let xv = linspace(-0.5,1,300), y = 1e-5im
    cal = [f1_II(x+1im*y) for x in xv]
    plot( xv, real.(cal))
    plot!(xv, imag.(cal))
    plot!(xv,  abs.(cal))
end

using Fitting
# The Pole Position
to_be_minimized(x,y) = abs(1./f1_II(x+1im*y))

let (v, pars, stat) = fit(to_be_minimized; starting_pars=[0.77^2, -0.15*0.77])
    sqrtsp = sqrt(complex(pars...))
    println("pole mass is $(real(sqrtsp)), pole width is $(-2*imag(sqrtsp))")
end
# plotlyjs()
# let xv = linspace(-1,1,100), yv = linspace(-0.5,0.5,100)
#     plot(layout=grid(2,1))
#     cal = [f1_I(x+1im*y) for y in yv, x in xv]
#     surfase!(xv, yv, log.(abs.(cal)),subplot=2)
#     cal = [f1_II(x+1im*y) for y in yv, x in xv]
#     plot!(xv, yv, log.(abs.(cal)), levels=40,subplot=1)
# end

let xv = linspace(0.1,2,100)
    plot(layout=grid(2,1))
    cal = [ f1_I(x) for x in xv]
    plot!(xv, real.(cal), subplot=1)
    plot!(xv, imag.(cal), subplot=1)
    cal = [ f1_I(x+1e-5im) for x in xv]
    plot!(xv, real.(cal), subplot=1)
    plot!(xv, imag.(cal), subplot=1)
    cal = [f1_II(x) for x in xv]
    plot!(xv, real.(cal), subplot=2)
    plot!(xv, imag.(cal), subplot=2)
    cal = [f1_II(x-1e-5im) for x in xv]
    plot!(xv, real.(cal), subplot=2)
    plot!(xv, imag.(cal), subplot=2)
end
