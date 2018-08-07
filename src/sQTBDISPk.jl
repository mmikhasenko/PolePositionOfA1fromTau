module sQTBDISPk

using JPACModel: ρ_qtb, ρ_qtb_i, ρ_qtb_disp, ρ_qtb_disp_i

export t_sqtb_dispk, t_sqtb_dispk_II, t_sqtb_dispk_II_external
export intens_sqtb_dispk

using Masses: mτ

# # --------------------------------------------------------------

K(s,m,g) = g^2/(m^2-s)
K(s,m,g,h) = g^2/(s*(m^2-s)+h)
K(s,m,g,mp,gp) = g^2/(m^2-s) + gp/(mp-s)
# --------------------------------------------------------------
function t_sqtb_dispk(s,args...)
    Kf = K(s,args...)
    1.0/(1.0/Kf-0.5im*s*ρ_qtb_disp(s))
end
function t_sqtb_dispk(s::Real,args...)
    Kf = K(s,args...)
    1.0/(1.0/Kf-0.5im*s*ρ_qtb_disp_i(s))
end
function t_sqtb_dispk_II(s,args...)
    Kf = K(s,args...)
    1.0/(1.0/Kf-1im*s*(ρ_qtb_disp(s)/2.0+ρ_qtb(s)))
end
function t_sqtb_dispk_II_external(qtb_external,s,args...)
    Kf = K(s,args...)
    1.0/(1.0/Kf-1im*s*(ρ_qtb_disp(s)/2.0+qtb_external(s)))
end


# ---------------------------------------------------------------------------------------
function intens_sqtb_dispk(s,pars...)
    A = t_sqtb_dispk(s,pars[1:end-1]...)
    phsp = ρ_qtb_i(s)
    abs2(pars[end]*(s*A))*phsp*(mτ^2-s)^2/s*(1+2*s/mτ^2)
 end

end
