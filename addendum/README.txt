Appendix C analytic-continuation investigation

Main deliverables: report.tex and report.pdf.
Numerical scratchpad: continuation_checks.py; results.json.
Figures: figures/; source paper and extracted text: sources/.

Reproduce with Python 3 and packages in requirements.txt:
  python continuation_checks.py
  pdflatex -interaction=nonstopmode -halt-on-error report.tex
  pdflatex -interaction=nonstopmode -halt-on-error report.tex
Run LaTeX from addendum/.

The numerical code implements a local physical-germ continuation, not a
certified global continuation to the a1 pole. It describes the baseline source at e508486; see docs/ for the subsequent Julia modernization.
Full-model results omit the common isobar normalization ns1; relative
comparisons and the analytic diagnosis do not depend on that constant.

Repository baseline: e508486ffa9370358e0039354be6449e603fa128
Source paper: https://arxiv.org/abs/1810.00016v2
Investigation date: 2026-10-15
