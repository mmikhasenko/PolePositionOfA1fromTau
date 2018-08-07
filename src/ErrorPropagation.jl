module ErrorPropagation

using Optim

export propagate_errors_to_BW_parameters, get_ellipce_from_covariance_from_mg
export get_ellipce_from_covariance, get_ellipce_from_the_set

"""
    propagate_errors_to_BW_parameters(spectral_over_gsq, pars, cov)

    the `spectral_over_gsq` is the spectral function devided by g^2.

    D(s) = m^2 - s - i g^2 F(s)

    The function returns Breit-Wigner'
        mass, width, error to the mass, error to the width, correlation.

"""
function propagate_errors_to_BW_parameters(spectral_over_gsq, pars, cov)
    m = pars[1]; g = pars[2];
    msq = m^2; gsq = g^2;
    B(s) = real(spectral_over_gsq(s))
    A(s) = imag(spectral_over_gsq(s))

    optimres = Optim.optimize(s->(msq+gsq*A(s)-s)^2, 1, 3)
    s0 = Optim.minimizer(optimres)
    e0 = sqrt(s0);
    # values
    A0 = A(s0); B0 = B(s0);
    # derivatives
    Bprime = (B(s0+1e-3)-B(s0-1e-3))/2e-3
    Aprime = (A(s0+1e-3)-A(s0-1e-3))/2e-3
    # ds0 / dg and ds0 / dm
    ds!dm = 2m   /(1-gsq*Aprime)
    ds!dg = 2g*A0/(1-gsq*Aprime)
    de!dm = ds!dm/(2*e0)
    de!dg = ds!dg/(2*e0)
    # dΓ!dg and dΓ!dm
    dΓ!dm = (-B0+m*Bprime*ds!dm)*gsq/msq
    dΓ!dg = (2B0+g*Bprime*ds!dg)*g/m
    Γ = gsq*B0/e0
    # finally combine
    de = sqrt(de!dm^2*cov[1,1]+de!dg^2*cov[2,2]+2*de!dg*de!dm*cov[1,2])
    dΓ = sqrt(dΓ!dm^2*cov[1,1]+dΓ!dg^2*cov[2,2]+2*dΓ!dg*dΓ!dm*cov[1,2])
    cv = dΓ!dm*de!dm*cov[1,1] + dΓ!dg*de!dg*cov[2,2] + (dΓ!dg*de!dm+dΓ!dm*de!dg)*cov[1,2]
    return e0, Γ, de, dΓ, cv/(de*dΓ)
#     return ds!dm, ds!dg
end

"""
    get_ellipce_from_covariance(cov,mn)

    bases on the covariance matrix cov it calculates eigenvalues and
    construct 2σ-ellipces around mean values given by mn.

"""
function get_ellipce_from_covariance(cov,mn)
#     cov_mΓ,mean_mΓ
    # eigen velue and eigenvectors
    eigval_mΓ = eigvals(cov);
    eigvec_mΓ = eigvecs(cov);
    # create an ellipse in the sky
    ϕ=linspace(0,2π,1000);
    x = sin.(ϕ); y = cos.(ϕ);
    chi2v_3sigma = 9.21
    chi2v_2sigma = 5.99
    x *= sqrt(eigval_mΓ[1]*chi2v_2sigma)
    y *= sqrt(eigval_mΓ[2]*chi2v_2sigma)
    θ = atan2(eigvec_mΓ[1,2],eigvec_mΓ[1,1])
    x, y = x*cos(θ)-y*sin(θ),x*sin(θ)+y*cos(θ)
    x, y = x+mn[1], y+mn[2]
    collect(zip(x, y))
end

"""
    get_ellipce_axis_from_covariance(cov,mn)

    bases on the covariance matrix cov it calculates eigenvalues and
    construct 2σ-ellipces around mean values given by mn.

"""
function get_ellipce_axis_from_covariance(cov,mn)
#     cov_mΓ,mean_mΓ
    # eigen velue and eigenvectors
    eigval_mΓ = eigvals(cov);
    eigvec_mΓ = eigvecs(cov);
    # create an ellipse in the sky
    x = [-1,0,0,1];
    y = [0,1,-1,0];
    chi2v_3sigma = 9.21
    chi2v_2sigma = 5.99
    x *= sqrt(eigval_mΓ[1]*chi2v_2sigma)
    y *= sqrt(eigval_mΓ[2]*chi2v_2sigma)
    θ = atan2(eigvec_mΓ[1,2],eigvec_mΓ[1,1])
    x, y = x*cos(θ)-y*sin(θ),x*sin(θ)+y*cos(θ)
    x, y = x+mn[1], y+mn[2]
    collect(zip(x, y))
end


"""
    get_ellipce_from_covariance_from_mg(spectral_over_gsq, pars, cov)

    simply combines the functions
     * `propagate_errors_to_BW_parameters(spectral_over_gsq, pars, cov)`
     * `get_ellipce_from_covariance(cov,mn)`

"""
function get_ellipce_from_covariance_from_mg(spectral_over_gsq, pars, cov)
    (e, Γ, de, dΓ, cor) = propagate_errors_to_BW_parameters(spectral_over_gsq, pars, cov)
    get_ellipce_from_covariance([de^2 cor*de*dΓ;cor*de*dΓ dΓ^2], [e, Γ])
end



"""
    get_ellipce_from_the_set(A,B)

    Construct covariance ellipces by the given arrays of `A` and `B`.
    It calculates covariance matrix and then calls
    `get_ellipce_from_covariance(cov,mn)`

"""
function get_ellipce_from_the_set(marr,Γarr)
    mΓ_mat = hcat(marr,Γarr)
    cov_mΓ = [cov(mΓ_mat[:,i],mΓ_mat[:,j]) for i=1:2,j=1:2]
    mean_mΓ = [mean(mΓ_mat[:,i]) for i=1:2]
    get_ellipce_from_covariance(cov_mΓ, mean_mΓ)
end

end
