module JPACModel

# registered packages
using Cuba
# my packages
push!(LOAD_PATH, "src")
using Isobars: f1_I, f1_II
using Masses: mπ, mπ2, mρ, mτ
using DalitzPlotAnalysis
using SomeRoutines

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
    s3_p, s3_m = 2mπ2+EE4/2+sqrtλλ/(2s1)*[1.0,-1.0]
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
    s3_p, s3_m = 2mπ2+EE4/2+sqrtλλ/(2s1)*[1.0,-1.0]
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
  [integrand([x,y]) for y in linspace(0,1,N), x in linspace(0,1,N)]
end

symm(args...) = qtb(args...) + interference(args...)
ρ_symm = symm
ρ_qtb = qtb

# ----------------------------------------------------------------------
# Dispersive integral
for name in ["symm", "qtb"]
    _ρ_name         = Symbol("ρ_"*name);
    _ρ_name_i       = Symbol("ρ_"*name*"_i");
    _ρ_name_disp0   = Symbol("ρ_"*name*"_disp0");
    _ρ_name_disp_r0 = Symbol("ρ_"*name*"_disp_r0");
    _ρ_name_disp    = Symbol("ρ_"*name*"_disp");
    _ρ_name_disp_i  = Symbol("ρ_"*name*"_disp_i");

    # # # # First lookup table for the ρsymm
    @eval $(_ρ_name_i) = make_lookup_with_asymptotics_and_tan_map($(_ρ_name), 9mπ2, 1e4)

    # # # Second, the dispersive integral
    # real argument
    @eval $(_ρ_name_disp0)(s::Real) = s/(π*1im)*quadgk(t->
        begin
            sp = 9mπ2+tan(t)
            $(_ρ_name_i)(sp)/(sp*(sp-s-1e-7im)*cos(t)^2)
        end
        ,0,π/2)[1]
    #shift
    @eval const $(_ρ_name_disp_r0) = imag($(_ρ_name_disp0)((mπ+mρ)^2))
    @eval $(_ρ_name_disp)(s::Real) = $(_ρ_name_disp0)(s)-1im*$(_ρ_name_disp_r0)

    # complex argument
    @eval function $(_ρ_name_disp)(s)
        (imag(s) == 0.0) && return $(_ρ_name_disp)(real(s))
        integr = quadgk(t->
            begin
                sp = 9mπ2+tan(t)
                $(_ρ_name_i)(sp)/(sp*(sp-s)*cos(t)^2)
            end
            ,0,π/2)[1]
        -1im*$(_ρ_name_disp_r0)+s/(π*1im)*integr
    end

    @eval $(_ρ_name_disp_i) = make_lookup_function($(_ρ_name_disp), linspace(9mπ2,3.5,60))
end

using Plots

plot(x->ρ_qtb(x), linspace(9mπ2, 3, 100))
plot( x->real(ρ_qtb_disp(x)), linspace(9mπ2, 3, 100))
plot!(x->imag(ρ_qtb_disp(x)), linspace(9mπ2, 3, 100))

# ----------------------------------------------------------------------
# ρ_qtb_i = make_lookup_with_asymptotics_and_tan_map(qtb, 9mπ2, 1e4)

# # # # Second, the dispersive integral
# # real argument
# ρ_qtb_disp0(s::Real)=s/π*quadgk(t->
#     begin
#         sp = 9mπ2+tan(t)
#         ρ_qtb_i(sp)/(sp*(sp-s-1e-7im)*cos(t)^2)
#     end
#     ,0,π/2)[1]
# const ρqr0 = real(ρtilde0((mπ+mρ)^2))
# ρ_qtb_disp(s::Real) = ρ_qtb_disp0(s)-ρqr0
#
# # complex argument
# ρ_qtb_disp(s::Complex{Float64})=-ρqr0+s/π*quadgk(t->
#     begin
#         sp = 9mπ2+tan(t)
#         ρi(sp)/(sp*(sp-s)*cos(t)^2)
#     end
#     ,0,π/2)[1]
#
# ρ_qtb_disp_i = make_lookup_function(ρ_qtb_disp, linspace(9mπ2,3.5,60))

# ------------------------------------------------------------------------------------
# amplitudes for the complex argument
for name in ["symm", "symm_disp", "qtb_disp", "qtb"]
    _t_name = Symbol("t_"*name)
    _ρ_name = Symbol("ρ_"*name)
    _ρ_name_i = Symbol("ρ_"*name*"_i")
    @eval function $(_t_name)(s,m,g)
        g^2/(m^2-s-0.5im*g^2*$(_ρ_name)(s))
    end
    @eval function $(_t_name)(s::Real,m,g)
        g^2/(m^2-s-0.5im*g^2*$(_ρ_name_i)(s))
    end
end

t_qtb(1.1, 1.1, 7)


# function t_symm(s::Real,m,g)
#     1.0/(m^2-s-0.5im*g^2*ρsi(s))
# end
# function t_symm(s,m,g)
#     1.0/(m^2-s-0.5im*g^2*symm(s))
# end


# ------------------------------------------------------------------------------------
# function t_symm_disp(s,m,g)
#     1.0/(m^2-s-g^2*ρstilde(s)/2.0)
# end
#
# function t_symm_disp(s::Real,m,g)
#     1.0/(m^2-s-g^2*ρstildei(s)/2.0)
# end

function t_symm_disp_II(s,m,g)
    g^2/(m^2-s-1im*g^2*(ρ_symm_disp(s)/2.0+symm(s)))
end

function t_symm_disp_II_external(s,m,g,symm_external)
    g^2/(m^2-s-1im*g^2*(ρ_symm_disp(s)/2.0+symm_external(s)))
end

# ------------------------------------------------------------------------------------
# function t_qtb_disp(s,m,g)
#     1.0/(m^2-s-g^2*ρtilde(s)/2.0)
# end
#
# function t_qtb_disp(s::Real,m,g)
#     1.0/(m^2-s-g^2*ρtildei(s)/2.0)
# end

function t_qtb_disp_II(s,m,g)
    g^2/(m^2-s-1im*g^2*(ρ_qtb_disp(s)/2.0+ρ_qtb(s)))
end

function t_qtb_disp_II_external(s,m,g,qtb_external)
    g^2/(m^2-s-1im*g^2*(ρ_qtb_disp(s)/2.0+qtb_external(s)))
end
# # ------------------------------------------------------------------------------------
# function t_qtb_disp_sqrts(s,m,g)
#     1.0/(m^2-s-g^2*ρtilde(s)/2.0)
# end
#
# function t_qtb_disp_sqrts(s::Real,m,g)
#     1.0/(m^2-s-g^2*ρtildei(s)/2.0)
# end
#
# function t_qtb_disp_sqrts_II(s,m,g)
#     1.0/(m^2-s-g^2*(ρtilde(s)/2.0+1.0im*qtb(s)))
# end

# ------------------------------------------------------------------------------------
# for the following functions the model and data must be defined
for model in ["symm", "symm_disp"]
      @eval function $(Symbol("intens_"*model))(s,pars...)
            A = $(Symbol("t_"*model))(s,pars[1:end-1]...)
            phsp = ρ_symm_i(s)
            abs2(pars[end]*A)*phsp*(mτ^2-s)^2/s*(1+2*s/mτ^2)
        end
end

for model in ["qtb", "qtb_disp"]
      @eval function $(Symbol("intens_"*model))(s,pars...)
            A = $(Symbol("t_"*model))(s,pars[1:end-1]...)
            phsp = ρ_qtb_i(s)
            abs2(pars[end]*A)*phsp*(mτ^2-s)^2/s*(1+2*s/mτ^2)
        end
end

let
    plot( x->real(t_qtb(x, 1.3, 7.0)), 0.3:0.03:3)
    plot!(x->real(t_symm(x, 1.3, 7.0)), 0.3:0.03:3)
    plot!( x->real(t_qtb_disp(x, 1.3, 7.0)), 0.3:0.03:3)
    plot!(x->real(t_symm_disp(x, 1.3, 7.0)), 0.3:0.03:3)

    plot!( x->imag(t_qtb(x, 1.3, 7.0)), 0.3:0.03:3)
    plot!(x->imag(t_symm(x, 1.3, 7.0)), 0.3:0.03:3)
    plot!( x->imag(t_qtb_disp(x, 1.3, 7.0)), 0.3:0.03:3)
    plot!(x->imag(t_symm_disp(x, 1.3, 7.0)), 0.3:0.03:3)
end

intens_qtb(1.1, 1.1, 7, 1.0)
intens_qtb_disp(1.1, 1.1, 7, 1.0)
intens_symm(1.1, 1.1, 7, 1.0)
intens_symm_disp(1.1, 1.1, 7, 1.0)

end
