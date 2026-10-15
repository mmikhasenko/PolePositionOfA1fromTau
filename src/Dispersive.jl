module Dispersive
using QuadGK: quadgk
export DispersiveLineshape, selfenergy, amplitude, intensity

"""
    DispersiveLineshape(rho; threshold, subtraction, rtol=1e-7)

Once-subtracted dispersive model for a real spectral density `rho(x)` supported
above `threshold`. All energy arguments are squared energies in GeV².
The real part of the self-energy vanishes at `subtraction`. The density must
be continuous, piecewise differentiable, and grow slower than x at infinity.
See the generated equations page for normalization and sheet conventions.
"""
struct DispersiveLineshape{F}
    rho::F
    threshold::Float64
    subtraction::Float64
    offset::Float64
    rtol::Float64
end

# The subtracted integrand has a removable singularity at x=s.
function raw_boundary(rho, a, s, rtol)
    s == 0 && return 0.0im
    if s < a
        return complex(s/π*quadgk(x->rho(x)/(x*(x-s)), a, Inf; rtol)[1])
    end
    if s == a
        rho(a) == 0 || throw(DomainError(s, "Nonzero density at threshold has a logarithmic singularity"))
        return complex(s/π*quadgk(x->rho(x)/(x*(x-s)), a, Inf; rtol)[1])
    end
    r = rho(s)
    function regular(x)
        if abs(x-s) < 1e-7*max(1,abs(s))
            h = min(1e-5*max(1,abs(s)), (s-a)/4)
            return (rho(s+h)-rho(s-h))/(2h)
        end
        s*(rho(x)-r)/(x*(x-s))
    end
    pv = quadgk(regular, a, s, Inf; rtol, atol=rtol*1e-3)[1]
    complex((pv+r*log(a/abs(a-s)))/π, r)
end
function DispersiveLineshape(rho; threshold, subtraction, rtol=1e-7)
    threshold > 0 || throw(ArgumentError("threshold must be positive"))
    rtol > 0 || throw(ArgumentError("rtol must be positive"))
    offset = real(raw_boundary(rho, threshold, subtraction, rtol))
    DispersiveLineshape(rho, Float64(threshold), Float64(subtraction), offset, Float64(rtol))
end

"""
    selfenergy(model, s; sheet=:I, continued_density=nothing)

Return Π(s), with Im Π(s+i0)=rho(s). A real argument means the upper rim.
For `sheet=:II`, supply a holomorphic `continued_density(s)` with the intended
contour history; ΠII=ΠI+2i*rho_cont is the adjacent-sheet formula below the cut.
No principal-root continuation of a physical squared amplitude is inferred.
"""
function selfenergy(model::DispersiveLineshape, s; sheet=:I, continued_density=nothing)
    a = model.threshold
    value = if isreal(s)
        raw_boundary(model.rho, a, real(s), model.rtol)-model.offset
    else
        s/π*quadgk(x->model.rho(x)/(x*(x-s)), a, Inf; rtol=model.rtol)[1]-model.offset
    end
    sheet == :I && return value
    sheet == :II || throw(ArgumentError("sheet must be :I or :II"))
    imag(s) < 0 || throw(ArgumentError("sheet II API expects Im(s)<0"))
    isnothing(continued_density) && throw(ArgumentError("sheet II requires continued_density"))
    value + 2im*continued_density(s)
end

"""
    amplitude(model, s, m, g; kwargs...)

Return g²/(m²-s-g²Π(s)/2), with m and g in GeV for a dimensionless density.
Keyword arguments are passed to `selfenergy`.
"""
amplitude(model::DispersiveLineshape, s, m, g; kwargs...) =
    g^2/(m^2-s-g^2*selfenergy(model,s; kwargs...)/2)

"""
    intensity(model, s, m, g, normalization=1; mtau=1.776)

Tau-decay spectrum |normalization*t(s)|² rho(s) (mtau²-s)²/s (1+2s/mtau²).
Returns zero outside the physical interval. Overall units depend on normalization.
"""
function intensity(model::DispersiveLineshape, s::Real, m, g, normalization=1; mtau=1.776)
    model.threshold < s < mtau^2 || return 0.0
    abs2(normalization*amplitude(model,s,m,g))*model.rho(s)*(mtau^2-s)^2/s*(1+2s/mtau^2)
end
end
