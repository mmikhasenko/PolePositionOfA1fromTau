using QuadGK
using Plots
theme(:default)
# plot(ρ,0,3)
# plot(x->real(BW(x,1.2,0.1)),0,3)
# plot!(x->imag(BW(x,1.2,0.1)),0,3)

function constructBW(mξ, Γξ, mi, gi)

    function U(σ)
        m = mξ
        Γ = Γξ # 0.15
        1.0/((m^2-σ)^2 + (m*Γ)^2)*sqrt(λ(σ,0,0))/σ
    end

    λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
    function ρ0(s)
        xf = s
        xm = 10*real(s)+imag(s)/100.0im
        quadgk(σ->U(σ)*sqrt(s-σ)/sqrt(s), 0.1, xm, xf)[1]
    end

    const ρn = ρ0(100.0)
    ρ(s) = ρ0(s)/ρn/(8π)

    function BW(s)
        1.0/(mi^2-s-1.0im*gi^2*ρ(s))
    end

    return BW
end

function constructBW(mξ, mi, gi)

    function BW(s)
        mπ = 0.14
        1.0/(mi^2-s-1.0im*gi^2*
            sqrt(1im*(s-(mξ+mπ)^2))*sqrt(s-(mξ-mπ)^2)*cis(-π/4)/
            (8π*s))
    end

    return BW
end

function plot_poles(amp)
    my_amp = amp
    @time let sxv = linspace(1e-3, 2.0, 100), syv = linspace(-1.0, 0.1, 100)
        cal = [my_amp(sx+1.0im*sy) for sy in syv, sx in sxv]
        contour(sxv,syv, log.(abs.(cal)), levels=100)
    end
end

# constuctBW(0.77,0.12,1.1,5.1)

plot_poles(constructBW(0.77,0.12,1.1,5.1))
plot_poles(constructBW(0.77,1.1,5))

using InteractNext

Dict("sin"=>:sin, "cos"=>:cos)


@manipulate for x in 1:0.001:3
    f = s->sin(s*x)
    plot(f,0,2)
end

#
@manipulate for g in linspace(1,6,10), Γξ in linspace(0.01,0.15,3), mξ in linspace(0.15,0.9,10), m in linspace(0.5,1.5,10)
    plot_poles(constructBW(mξ,Γξ,m,g))
end

@manipulate for g in linspace(1,6,10), mξ in linspace(0.15,0.9,10), m in linspace(0.5,1.5,10)
    plot_poles(constructBW(mξ,m,g))
end
