module JPACModel

# registered packages
using Cuba
# my packages
using QuadGK: quadgk
using ..Dispersive: DispersiveLineshape, selfenergy
using ..Isobars: f1_I, f1_II
using ..Masses: mπ, mπ2, mρ, mτ
using ..DalitzPlotAnalysis
using ..SomeRoutines

export qtb, interference, symm
export interference_array

for name in ["symm", "symm_disp", "qtb_disp", "qtb"]
    @eval export $(Symbol("ρ_"*name))
    @eval export $(Symbol("ρ_"*name*"_i"))
    @eval export $(Symbol("t_"*name))
    @eval export $(Symbol("intens_"*name))
end

export t_symm_disp_II, t_symm_disp_II_external
export t_qtb_disp_II, t_qtb_disp_II_external
# -----

# The special implementation for the Float64
# Suppose to work way faster

# function qtb(s::Float64)
#     (s ≤ 9mπ2) && return 0.0
#
#     function integrand(x,f)
#         s1min = 4.0 * mπ2
#         s1max = (sqrt(s) - mπ)^2
#
#         s1 =  s1min + x[1] * (s1max - s1min)
#         f[1] = abs2(f1_I(s1)) * sqrt(λ(s,s1,mπ2)*λ(s1,mπ2,mπ2)) / (s*s1)*(s1max-s1min)
#     end
#     return vegas(integrand, 1, 1)[1][1]/(8π)^2/(2π)
# end

# function qtb(s)
#   function integrand(x,f)
#     s1min = 4.0 * mπ2
#     s1max = (sqrt(s) - mπ)^2
#
#     s1 =  s1min + x[1] * (s1max - s1min)
#     f[1], f[2] = reim(f1_I(s1)*f1_II(s1) * sqrt(λ(s,s1,mπ2)*λ(s1,mπ2,mπ2)) / (s*s1)*(s1max-s1min))
#   end
#   return complex(vegas(integrand, 1, 2)[1]...)/(8π)^2/(2π)
# end

function qtb(s::Real)
  s <= 9mπ2 && return 0.0
  s1min = 4.0 * mπ2
  s1max = (sqrt(s) - mπ)^2
  function integrand(s1)
    λλ = λ(s,s1,mπ2)*λ(s1,mπ2,mπ2)
    (λλ < 0.0) && return 0.0
    abs2(f1_I(s1)) * sqrt(λλ) / (s*s1)
  end
  return quadgk(integrand,s1min,s1max)[1]/(8π)^2/(2π)
end

function qtb(s, rr = 1.0)
  s1min = 4.0 * mπ2
  sqrts = sqrt(s)
  s1max = (sqrts - mπ)^2
  s1mid = (rr > 0.0) ? rr*real(s1max)+1e-4im*sign(imag(s)) : s1min
  return quadgk(s1->f1_I(s1)*f1_II(s1) *
                     sqrt((sqrts+mπ)^2-s1) * sqrt((sqrts-mπ)^2-s1) *
                     sqrt(s1-4mπ2) * sqrt(s1) / (s*s1),
    s1min,s1mid,s1max)[1]/(8π)^2/(2π)
end


# qtb(1.1)
# ------------------------------------------------------------

function angular_pow2(s, s1, z, m1sq, m2sq, m3sq)
    θ23 = acos(z)
    θ3 = acos(cross_basis_cosθ3(s1, z, m1sq, m2sq, m3sq, s))
    θ12 = acos(cross_basis_cosθ12(s1, z, m1sq, m2sq, m3sq, s))
    # the minus sign is due to the (-1)^S for isospin
    res = -Wignerd(1,0,0,θ3+θ12-θ23);
    return res
end

function angular_pow2(s, s1, s3)
    sqrts = sqrt(s); sqrts1 = sqrt(s1); sqrts3 = sqrt(s3);
    topW = s1*s3*(2*s*(mπ2 + s - s1 - s3) -
        (mπ2 + s - s1)*(mπ2 + s - s3))*(3*mπ2 + s - 2*s1 - s3)*(3*mπ2 + s - s1 - 2*s3) -
        4*(sqrts*sqrts1*s3*(3*mπ2 + s - 2*s1 - s3) +
        sqrts*s1*sqrts3*(3*mπ2 + s - s1 - 2*s3) +
        sqrts1*sqrts3*(-2*s*(mπ2 + s - s1 - s3) +
        (mπ2 + s - s1)*(mπ2 + s - s3))) *
            (mπ2^3 - 2*mπ2^2*s + mπ2*s^2 - 3*mπ2*s1*s3 - s*s1*s3 + s1^2*s3 + s1*s3^2)
    bottomW = sqrts1*sqrt(s1-4mπ2)*sqrts3*sqrt(s3-4mπ2)*λ(s,s1,mπ2)*λ(s,s3,mπ2)
    return -topW/bottomW
end

function angular_pow2(s, s1, z, m1sq, m2sq, m3sq, s3)
    θ23 = acos(z)
    cosθ3 = cross_basis_cosθ3(s1, z, m1sq, m2sq, m3sq, s, s3)
    θ3 = acos(cosθ3) # this causes problem
    cosθ12 = cross_basis_cosθ12(s1, z, m1sq, m2sq, m3sq, s, s3)
    θ12 = acos(cosθ12)
    # the minus sign is due to the (-1)^S for isospin
    # return z
    # return acos(cosθ3)
    # return acos(cosθ12) #z + acos(cosθ3) +
    res = -Wignerd(1,0,0,θ3+θ12-θ23);
    return res
end


function interference(s::Float64)
    (s ≤ 9mπ2) && return 0.0
    function integrand(x,f)
      s1min = 4.0 * mπ2
      s1max = (sqrt(s) - mπ)^2

      s1 =  s1min + x[1] * (s1max - s1min)
      z = -1.0 + 2.0*x[2]
      s3 = cross_basis_s3(s1, z, mπ2, mπ2, mπ2, s)
      ampl2 = real(f1_I(s1)*f1_II(s3))*angular_pow2(s, s1, s3) #z, mπ2, mπ2, mπ2
      f[1] = ampl2 * sqrt(λ(s,s1,mπ2)*λ(s1,mπ2,mπ2)) / (s*s1)*(s1max-s1min)
    end
    return (cuhre(integrand, 2, 1)[1][1]/(8π)^2/(2π))::Float64
end

function interference(s, rr = 1.0)
  sqrts = sqrt(s)
  function integrand(x,f)
    s1min = 4mπ2
    s1max = (sqrts - mπ)^2

    s1med = (rr > 0.0) ? rr*real(s1max)+1e-3im*imag(s1max) : s1max
    s1 = zero(s)
    ampl2 = one(s)
    if x[1] < 0.5
        s1 = s1min + x[1] * 2*(s1med - s1min);
        ampl2 =             2*(s1med - s1min)
    else
        s1 = s1med + (x[1]-0.5) * 2*(s1max - s1med);
        ampl2 =                   2*(s1max - s1med);
    end
        # s3
    EE4 = s-mπ2-s1
    sqrtλλ = sqrt((sqrts-mπ)^2-s1)*sqrt((sqrts+mπ)^2-s1)*sqrt(s1-4mπ2)*sqrt(s1)
    s3_p, s3_m = (2mπ2+EE4/2+sqrtλλ/(2s1), 2mπ2+EE4/2-sqrtλλ/(2s1))
    s3_c =5mπ2#real(s3_p+s3_m)/2

    s3 = zero(s)
    if x[2] < 0.5
        s3 = s3_m + x[2] * 2*(s3_c - s3_m);
        ampl2 *=           2*(s3_c - s3_m)
    else
        s3 = s3_c + (x[2]-0.5) * 2*(s3_p - s3_c);
        ampl2 *=                 2*(s3_p - s3_c);
    end

    if (s1 ≈ 4mπ2) || (s1 ≈ (sqrts-mπ)^2) || (s3 ≈ 4mπ2) || (s3 ≈ (sqrts-mπ)^2)
      f[1],f[2] = 0.0, 0.0
      return
    end

    ampl2 *= f1_I(s1)*f1_II(s3)*angular_pow2(s, s1, s3)

    f[1],f[2] = reim(ampl2)
  end
  return complex(cuhre(integrand, 2, 2)[1]...)/(8π)^2/(2π*s) #, abstol=1e-40, reltol=1e-42
end

function interference_array(s, rr = 1.0, N=10)
  sqrts = sqrt(s)
  function integrand(x)
    s1min = 4mπ2
    s1max = (sqrts - mπ)^2

    s1med = (rr > 0.0) ? rr*real(s1max)+1e-3im*imag(s1max) : s1max
    s1 = zero(s)
    ampl2 = one(s)
    if x[1] < 0.5
        s1 = s1min + x[1] * 2*(s1med - s1min);
        ampl2 =             2*(s1med - s1min)
    else
        s1 = s1med + (x[1]-0.5) * 2*(s1max - s1med);
        ampl2 =                   2*(s1max - s1med);
    end
        # s3
    EE4 = s-mπ2-s1
    sqrtλλ = sqrt((sqrts-mπ)^2-s1)*sqrt((sqrts+mπ)^2-s1)*sqrt(s1-4mπ2)*sqrt(s1)
    s3_p, s3_m = (2mπ2+EE4/2+sqrtλλ/(2s1), 2mπ2+EE4/2-sqrtλλ/(2s1))
    s3_c = 6mπ2#real(s3_p+s3_m)/2

    s3 = zero(s)
    if x[2] < 0.5
        s3 = s3_m + x[2] * 2*(s3_c - s3_m);
        ampl2 *=           2*(s3_c - s3_m)
    else
        s3 = s3_c + (x[2]-0.5) * 2*(s3_p - s3_c);
        ampl2 *=                 2*(s3_p - s3_c);
    end

    if (s1 ≈ 4mπ2) || (s1 ≈ (sqrts-mπ)^2) || (s3 ≈ 4mπ2) || (s3 ≈ (sqrts-mπ)^2)
      return 0.0im
    end

    ampl2 *= f1_I(s1)*f1_II(s3) * angular_pow2(s, s1, s3)

    return ampl2
  end
  [integrand([x,y]) for y in range(0,1,N), x in range(0,1,N)]
end

symm(args...) = qtb(args...) + interference(args...)
ρ_symm = symm
ρ_qtb = qtb


# Build expensive interpolation tables only when requested, never at import time.
const density_cache = Dict{Symbol,Any}()
const model_cache = Dict{Symbol,Any}()
function density(name::Symbol)
    get!(density_cache, name) do
        f = name == :qtb ? qtb : name == :symm ? symm : error("Unknown density: $name")
        make_lookup_with_asymptotics_and_tan_map(f, 9mπ2, 1e4)
    end
end
function dispersive_model(name::Symbol)
    get!(model_cache, name) do
        DispersiveLineshape(density(name); threshold=9mπ2, subtraction=(mπ+mρ)^2)
    end
end
for name in (:qtb, :symm)
    rhoi = Symbol("ρ_", name, "_i")
    rhod = Symbol("ρ_", name, "_disp")
    rhodi = Symbol("ρ_", name, "_disp_i")
    @eval $rhoi(s::Real) = density($(QuoteNode(name)))(s)
    @eval $rhod(s) = -im*selfenergy(dispersive_model($(QuoteNode(name))), s)
    @eval $rhodi(s::Real) = $rhod(s)
end
for name in (:symm, :symm_disp, :qtb_disp, :qtb)
    t = Symbol("t_", name)
    rho = Symbol("ρ_", name)
    rhoi = Symbol("ρ_", name, "_i")
    @eval $t(s,m,g) = g^2/(m^2-s-0.5im*g^2*$rho(s))
    @eval $t(s::Real,m,g) = g^2/(m^2-s-0.5im*g^2*$rhoi(s))
    phase = startswith(string(name), "qtb") ? :ρ_qtb_i : :ρ_symm_i
    intensity = Symbol("intens_", name)
    @eval $intensity(s, m, g, n) = (9mπ2 < s < mτ^2) ?
        abs2(n*$t(s,m,g))*$phase(s)*(mτ^2-s)^2/s*(1+2s/mτ^2) : zero(m+g+n)
end
for name in (:qtb, :symm)
    tII = Symbol("t_", name, "_disp_II")
    external = Symbol(tII, "_external")
    rhod = Symbol("ρ_", name, "_disp")
    rho = Symbol("ρ_", name)
    @eval $external(s,m,g,continued) = g^2/(m^2-s-im*g^2*($rhod(s)/2+continued(s)))
    @eval $tII(s,m,g) = $external(s,m,g,$rho)
end
end
