module FunctionsRho

push!(LOAD_PATH, pwd()*"../modules")
# my modules
using isobars
using DalitzPlotAnalysis

export

# general function should not be defined
ρ(s) = 0
# real agrument
function ρ(s::Real)
    (s ≤ 9mπ2) && return 0.0
    1/(8*π*s)/(2π)*quadgk(s1->begin
        sqrt(λ(s,s1,mπ^2))*abs(f1_I(s1))^2*sqrt(λ(s1,mπ2,mπ2))/(8*π*s1)
        end, 4*mπ2,(sqrt(s)-mπ)^2)[1]
end

end
