using Documenter, PolePositionOfA1fromTau
include(joinpath(@__DIR__,"..","scripts","generate_docs.jl"))
generate_docs(check=true)
makedocs(
    sitename="Pole position of the a₁ from τ decay",
    modules=[PolePositionOfA1fromTau.Dispersive],
    format=Documenter.HTML(prettyurls=get(ENV,"CI","false")=="true", edit_link="main"),
    pages=["Home"=>"index.md", "Equations"=>"lineshape.md", "Scripts"=>"scripts.md",
           "API"=>"api.md", "Addendum"=>"addendum.md", "Citation"=>"citation.md"],
    doctest=true, checkdocs=:exports, warnonly=false,
)
