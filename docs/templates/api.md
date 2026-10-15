# API

```@docs
DispersiveLineshape
selfenergy
amplitude
intensity
qtb
interference
symm
```

The historical modules remain available below `PolePositionOfA1fromTau`: `Masses`, `Isobars`, `DalitzPlotAnalysis`, `JPACModel`, `SomeRoutines`, `Fitting`, `ErrorPropagation`, and `sQTBDISPk`.

Use explicit imports, for example `using PolePositionOfA1fromTau.JPACModel: t_qtb_disp`. The old practice of adding `src` to `LOAD_PATH` is no longer needed. Compatibility entry points under `tests/` delegate to the maintained scripts; automated tests live in `test/runtests.jl`.
