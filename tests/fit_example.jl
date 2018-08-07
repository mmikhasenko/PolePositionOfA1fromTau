push!(LOAD_PATH, "src")
using Fitting

data = hcat(collect(1:100),rand(100),[0.1 for i=1:100])
myf = build_chi2(((x,a,b)->a+b*x^2), data)

myf(1,1)

fit_with_grad(myf, starting_pars=[2,0.01])

fit(myf, starting_pars=[2,0.01])
