module Fitting

export build_chi2, fit, fit_with_grad

using NLopt
using ForwardDiff


"""
   build_chi2(f, data)

   The dataset `data` which is assumed to have three colums:
     - `x`-value, `y`-value, `dy`-error.

    If the function χ2 depends on n-arguments,
        the returned fuction will depend on `n-1`-parameters.

"""
function build_chi2(f, data)
    function mychi2(pars...)
        sum([(data[i,2]-f(data[i,1],pars...))^2/data[i,3]^2 for i in 1:size(data,1)])
    end
    return mychi2
end

"""
   build_chi2(f, data, cov))

   The dataset `data` which is assumed to have three colums:
     - `x`-value, `y`-value, `dy`-error.

   The covariance matrix `cov` should be a matrix of Ncol x Ncol size.

"""
function build_chi2(f, data, cov)
    function mychi2(pars...)
        Δ = data[:,2]-f.(data[:,1],pars...)
        inv_cov = inv(Symmetric(cov))
        transpose(Δ)*inv_cov*Δ
    end
    return mychi2
end

"""
   fit(χ2; starting_pars)

   The function minimizes `χ2` function starting from the point given by
       the `starting_pars`. It used NLopt library. The number and the order
       of parameters are deduced from the `starting_pars` vector.

"""
function fit(χ2; starting_pars::Vector{Float64} = error("Starting values are required!"),
                 algorithm::Symbol = :LN_COBYLA, verbose::Bool=true)
    to_minimize(x::Vector, grad::Vector) = let v = χ2(x...); verbose && @show v,x...; v end
    opt = Opt(algorithm, length(starting_pars))
    min_objective!(opt, to_minimize)
    xtol_rel!(opt,1e-6)
    optimize(opt, starting_pars)
end

"""
   fit_with_grad(χ2; starting_pars)

   The function minimizes `χ2` function starting from the point given by
   the `starting_pars`. The gradient is calculated automatically using ForwardDiff package.
   NLopt method 'LD_LBFGS' is called for the minimization. The number and the order
   of parameters are deduced from the `starting_pars` vector.

"""
function fit_with_grad(χ2; starting_pars::Vector{Float64} = error("Starting values are required!"),
                           algorithm::Symbol = :LD_MMA, verbose::Bool=true)
    function to_minimize(x::Vector, grad::Vector)
        v = χ2(x...);
        (length(grad) > 0) && (grad .= ForwardDiff.gradient(p->χ2(p...), x))
        verbose && (@show v,x...,grad)
        v
    end
    opt = Opt(algorithm, length(starting_pars)) # LD_LBFGS
    min_objective!(opt, to_minimize)
    xtol_rel!(opt,1e-6)
    (minf,pars,ret) = optimize(opt, starting_pars)

    # matrix of second derivatives
    val = ForwardDiff.hessian(p->χ2(p...),pars)
    cov_mat = inv(val./2.0)
    # return

    minf,pars,ret,cov_mat
end



end
