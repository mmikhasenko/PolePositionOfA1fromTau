# Analytic-continuation addendum

The former `research/` scratchpad is now [`addendum/`](https://github.com/mmikhasenko/PolePositionOfA1fromTau/tree/main/addendum).

Read the [PDF report](https://github.com/mmikhasenko/PolePositionOfA1fromTau/blob/main/addendum/report.pdf) or its [LaTeX source](https://github.com/mmikhasenko/PolePositionOfA1fromTau/blob/main/addendum/report.tex). It explains the original motivation for the second Appendix C contour, reproduces the endpoint plots, constructs the physical local continuation, checks limiting cases, and reviews related literature.

The local toy ambiguity is resolved by fixing the physical boundary data and tracking the two endpoint sectors. Global transport to the resonance pole remains unresolved: an undeformed contour can cross an isobar pole. The modern Julia package retains the historical complex-contour evaluator for investigation, but does not promote it to a certified continuation.

The report describes the baseline `e508486` implementation. Modernization replaces its eager tables, plotting on import, and fixed imaginary regulator; the contour discussion remains applicable.

To reproduce its calculations:

```sh
python3 -m venv /tmp/a1-addendum
/tmp/a1-addendum/bin/pip install -r addendum/requirements.txt
/tmp/a1-addendum/bin/python addendum/continuation_checks.py
/tmp/a1-addendum/bin/python addendum/pole_crossing_check.py
cd addendum
pdflatex -interaction=nonstopmode -halt-on-error report.tex
pdflatex -interaction=nonstopmode -halt-on-error report.tex
```

The source PDFs used by the figures are retained in `addendum/sources/`. Numerical experiments are independent Python checks; they are not a new ALEPH fit. Literature reading notes record which historical full texts could and could not be verified.
