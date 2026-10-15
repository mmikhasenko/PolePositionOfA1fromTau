module Dispersive
using QuadGK: quadgk
export DispersiveLineshape, selfenergy, amplitude, intensity

"""
    DispersiveLineshape(rho; threshold, subtraction, rtol=1e-7)

Construct a once-subtracted three-body self-energy model for use in a
Breit–Wigner-like amplitude. Its real part is obtained from a dispersion
integral over the supplied spectral density.

# Parameters
- `rho`: callable real, dimensionless spectral density `rho(x)`, supported
  above `threshold`; `x` is an invariant mass squared in GeV². For the paper's
  symmetrized a₁ → ρπ model, use the density prepared by `JPACModel.density(:symm)`.
- `threshold`: positive lower integration limit in GeV² (physically `9mπ²`).
- `subtraction`: squared energy in GeV² where `Re Π(subtraction+i0)=0`.
- `rtol`: positive relative quadrature tolerance; does not control density
  interpolation error.

The density must be continuous, piecewise differentiable, and grow slower
than `x` at infinity. Construction evaluates the real subtraction constant.
`JPACModel.dispersive_model(:symm)` prepares and caches the physical model.

# Example
```jldoctest
julia> using PolePositionOfA1fromTau

julia> model = DispersiveLineshape(x -> 1.0; threshold=0.2, subtraction=0.6);

julia> abs(real(selfenergy(model, 0.6))) < 1e-10
true
```
The constant density in this example is an analytic benchmark, not a₁ phase space.
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

Return the dimensionless subtracted loop Π(s), with `Im Π(s+i0)=model.rho(s)`.

# Parameters
- `model`: a [`DispersiveLineshape`](@ref), including the density and subtraction.
- `s`: real or complex three-pion invariant mass squared in GeV². A real
  argument denotes the upper rim of the physical cut.
- `sheet`: `:I` for the dispersive first sheet, or `:II` for its continuation
  into the adjacent lower-half-plane sheet (requires `imag(s)<0`).
- `continued_density`: callable holomorphic continuation of the density,
  required for `:II`; its contour and branch history must be supplied by the caller.

On sheet II the returned loop is `ΠI(s)+2im*continued_density(s)`.
The dispersion integral fixes the first sheet; it does not by itself certify
complex Dalitz-contour transport to a resonance pole.

# Example
```jldoctest
julia> using PolePositionOfA1fromTau

julia> model = PolePositionOfA1fromTau.JPACModel.dispersive_model(:symm);

julia> isapprox(imag(selfenergy(model, 1.2)), model.rho(1.2))
true
```
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
    amplitude(model, s, m, g; sheet=:I, continued_density=nothing)

Return the dimensionless Breit–Wigner-like amplitude
`t(s) = g²/(m²-s-g²Π(s)/2)`. For the symmetrized model, Π includes the
Dalitz-integrated a₁ → ρπ → 3π decay chains and their interference.

# Parameters
- `model`: a [`DispersiveLineshape`](@ref) specifying the self-energy.
- `s`: three-pion invariant mass squared in GeV²; real values mean the upper rim.
- `m`: bare resonance mass parameter in GeV, not the extracted pole mass.
- `g`: coupling in GeV for the dimensionless spectral-density convention.
- `sheet`, `continued_density`: passed to [`selfenergy`](@ref).

# Example
```jldoctest
julia> using PolePositionOfA1fromTau

julia> model = PolePositionOfA1fromTau.JPACModel.dispersive_model(:symm);

julia> t = amplitude(model, 1.2, 1.3, 7.0); # illustrative m and g

julia> isapprox(imag(inv(t)), -model.rho(1.2)/2)
true
```
"""
amplitude(model::DispersiveLineshape, s, m, g; kwargs...) =
    g^2/(m^2-s-g^2*selfenergy(model,s; kwargs...)/2)

"""
    intensity(model, s, m, g, normalization=1; mtau=1.776)

Return the real-axis tau-decay spectrum
`|normalization*t(s)|² rho(s) (mtau²-s)²/s (1+2s/mtau²)`.

# Parameters
- `model`: a [`DispersiveLineshape`](@ref); its density is also the final-state factor.
- `s`: real three-pion invariant mass squared in GeV².
- `m`, `g`: bare resonance mass and coupling in GeV, as in [`amplitude`](@ref).
- `normalization`: overall production-amplitude normalization (default `1`);
  its modulus squared scales the spectrum and determines its overall units.
- `mtau`: tau mass in GeV (default `1.776`, the historical model value).

Returns zero outside `model.threshold < s < mtau²`. This is a spectrum,
not a probability density normalized to unit integral.

# Example
```jldoctest
julia> using PolePositionOfA1fromTau

julia> model = PolePositionOfA1fromTau.JPACModel.dispersive_model(:symm);

julia> intensity(model, 1.2, 1.3, 7.0, 1.0) > 0 # illustrative parameters
true
```
"""
function intensity(model::DispersiveLineshape, s::Real, m, g, normalization=1; mtau=1.776)
    model.threshold < s < mtau^2 || return 0.0
    abs2(normalization*amplitude(model,s,m,g))*model.rho(s)*(mtau^2-s)^2/s*(1+2s/mtau^2)
end
end
