include("common.jl")
objective, data, indices = prepared_fit(:qtb)
samples = readdlm(joinpath(ROOT,"data","aleph_3pi_bstp.txt"))
# Bootstrap rows are the nonzero bins of the complete data table.
full = readdlm(joinpath(ROOT,"data","ALEPH.2005.3pi.formatted"))
nonzero = findall(!iszero,full[:,2])
positions = [findfirst(==(i),nonzero) for i in indices]
@assert size(samples,1) == length(nonzero)
results = [begin
    y = samples[positions,j]
    result = Fitting.fit(objective(y); starting_pars=[1.2,7.0,0.02],algorithm=:LN_NELDERMEAD,verbose=false)
    [j,result[1],result[2]...]
end for j in 1:min(size(samples,2),SMOKE ? 2 : 100)]
write_table("bootstrap_qtb.csv",["replica","chi2","m","g","normalization"],permutedims(hcat(results...)))
