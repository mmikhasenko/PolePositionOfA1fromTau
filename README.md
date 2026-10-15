# Pole position of the a₁ from τ decay

Julia 1.13 code for [Mikhasenko et al., Phys. Rev. D 98, 096021 (2018)](https://doi.org/10.1103/PhysRevD.98.096021). Copy the paper's BibTeX from [CITATION.bib](CITATION.bib).

```sh
julia +1.13 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
julia +1.13 --project=. scripts/export_lineshape.jl
julia +1.13 --project=. scripts/generate_docs.jl
julia +1.13 --project=. docs/make.jl
```

| Start here | Purpose |
|:--|:--|
| `scripts/export_lineshape.jl` | Dispersive lineshape CSV + TOML ([example](examples/lineshape_qtb.csv)) |
| `scripts/fit_ALEPH.jl` | Fit the bundled data with QTB/SYMM models |
| `scripts/ALEPH_bootstrap.jl` | Parameter replicas from the supplied bootstrap data |
| `scripts/export_phase_space.jl`, `scripts/plot_JPACModel.jl` | Phase-space tables and example plot |
| `scripts/precalculate_qtb_interference.jl` | Diagnostic complex grids |
| [Full script map](docs/src/scripts.md) | Remaining studies and migration notes |

[Equations and definitions](docs/src/lineshape.md) · [API](docs/src/api.md) · [Documentation](https://mmikhasenko.github.io/PolePositionOfA1fromTau/). Generated Markdown lives in `docs/src/`; edit `docs/templates/` and regenerate. Documenter builds HTML in `docs/build/`; CI checks tests, doctests, and generated pages, uploads its report, and publishes successful main builds to `gh-pages` (served by GitHub Pages from that branch). Script output goes to `artifacts/`; `--smoke` runs smaller examples.

The [analytic-continuation addendum](addendum/report.pdf) ([source](addendum/report.tex), formerly `research/`) explains Appendix C and reviews related papers. Its local resolution does **not** certify global contour transport to the pole. This modernization makes the real-axis model reusable; it does not claim a revised pole result.
