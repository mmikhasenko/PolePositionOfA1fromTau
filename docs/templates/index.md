# Pole position of the a₁ from τ decay

Julia 1.13 implementation accompanying [Mikhasenko et al., Physical Review D 98, 096021 (2018)](https://doi.org/10.1103/PhysRevD.98.096021).

## Reproduce locally

From the repository root:

```sh
julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia +1.13 --project=. scripts/export_lineshape.jl
julia +1.13 --project=. scripts/generate_docs.jl
julia +1.13 --project=. docs/make.jl
```

The committed `Manifest.toml` fixes the dependency graph. Markdown in `docs/src/` is generated from `docs/templates/`, the script inventory, and `CITATION.bib`; edit those inputs and regenerate. Use `scripts/generate_docs.jl --check` to detect stale pages. Documenter builds `docs/build/index.html`, checks references and doctests, and fails CI on errors. CI retains the HTML and build log as downloadable artifacts.

Read [the equations](lineshape.md), [the runnable scripts](scripts.md), [the API](api.md), and [the analytic-continuation addendum](addendum.md). CSV and TOML output goes to `artifacts/` (override with `A1_OUTPUT_DIR`). Example parameters are explicitly illustrative; the upgrade is not a refit of the published pole.

## Reusable model

```jldoctest
julia> using PolePositionOfA1fromTau

julia> model = DispersiveLineshape(x -> 1.0; threshold=0.2, subtraction=0.6);

julia> round(imag(selfenergy(model, 0.7)); digits=6)
1.0

julia> round(imag(inv(amplitude(model, 0.7, 1.3, 7.0))); digits=6)
-0.5
```

For the paper's QTB spectral density, use

```julia
using PolePositionOfA1fromTau: JPACModel
model = JPACModel.dispersive_model(:qtb) # :symm includes real-axis interference
value = amplitude(model, 1.2, 1.3, 7.0)
```

The first call builds a density lookup table. Importing the package does not fit data, create plots, or build these tables.
