# Lineshape and dispersive integral

The lineshape is a Breit–Wigner-like parametrization whose energy-dependent self-energy describes a₁ → ρπ → 3π, symmetrized over the three-pion decay Dalitz plot. The Dalitz integral supplies the absorptive part; a subtracted dispersion integral supplies the corresponding real part and implements the first-sheet analytic structure of the three-body loop. The QTB variant omits interference between the two ρπ decay chains.

These formulas specify the implemented convention without requiring the historical plotting scripts. The real-axis model follows the paper's QTB and symmetrized phase-space construction; complex interference remains subject to the [addendum](addendum.md).

## Quantities and units

| Symbol | Meaning | Value / units |
|:--|:--|:--|
| ``s`` | total three-pion invariant mass squared | GeV² |
| ``u,v`` | two-pion invariant masses squared | GeV² |
| ``m_\pi,m_\rho,m_\tau`` | fixed historical masses | 0.13957, 0.77526, 1.776 GeV |
| ``\Gamma_\rho`` | rho on-shell width | 0.1491 GeV |
| ``R`` | rho barrier radius | 5 GeV⁻¹ |
| ``a`` | three-pion threshold ``9m_\pi^2`` | GeV² |
| ``s_0`` | subtraction point ``(m_\rho+m_\pi)^2`` | GeV² |
| ``m,g`` | fitted bare resonance mass and coupling | GeV |
| ``n`` | production normalization | sets spectrum units |
| ``\rho(s)`` | QTB or symmetrized spectral density | dimensionless |
| ``\Pi(s)`` | subtracted three-body loop | dimensionless |
| ``t(s)`` | resonance amplitude | dimensionless |

These fixed masses reproduce the repository model; they are not a claim to current PDG precision. ``m`` is a fit parameter, not the complex pole mass. Define

```math
\lambda(x,y,z)=x^2+y^2+z^2-2xy-2yz-2zx.
```

## Normalized rho isobar

For real ``u\ge4m_\pi^2``, define

```math
p(u)=\tfrac12\sqrt{u-4m_\pi^2},\quad p_0=p(m_\rho^2),\quad
r_2(u)=\frac{p(u)}{4\pi\sqrt u},\quad
F(u)=\frac{p(u)^2}{p_0^2}\frac{R^{-2}+p_0^2}{R^{-2}+p(u)^2},
```

```math
w(u)=m_\rho\Gamma_\rho\frac{r_2(u)}{r_2(m_\rho^2)}F(u),\qquad
f_{\rm raw}^{\mp}(u)=\frac{\sqrt{F(u)}}{m_\rho^2-u\mp i w(u)}.
```

The upper sign (minus in the denominator) is `Isobars.f1r_I`; the other factor is `f1r_II`. Set

```math
N=\int_{4m_\pi^2}^{\infty}\frac{du}{2\pi}\,
 f_{\rm raw}^{-}(u) f_{\rm raw}^{+}(u) r_2(u),\qquad
 f^{\mp}(u)=f_{\rm raw}^{\mp}(u)/\sqrt N.
```

Thus ``\int du\,f^-f^+r_2/(2\pi)=1``. ``f`` has units GeV⁻¹; ``N`` has units GeV⁻². On the physical axis ``f^+=\overline{f^-}``. Away from it, continue the two factors separately, including all roots; complex conjugation of a value at complex ``u`` is not an analytic operation.

## Spectral density

Let ``b=(\sqrt s-m_\pi)^2``. The quasi-two-body density is zero for ``s\le a`` and otherwise

```math
\rho_{\rm QTB}(s)=\frac{1}{2\pi(8\pi)^2}
\int_{4m_\pi^2}^{b}du\,|f^-(u)|^2
\frac{\sqrt{\lambda(s,u,m_\pi^2)\lambda(u,m_\pi^2,m_\pi^2)}}{su}.
```

For the symmetrized model define ``M=m_\pi^2``, ``x=\sqrt s``, ``y=\sqrt u``, ``z=\sqrt v``, and

```math
v_\pm(s,u)=\frac{s+3M-u}{2}\pm
 \frac{\sqrt{\lambda(s,u,M)\lambda(u,M,M)}}{2u}.
```

The following polynomial specifies the full implemented angular kernel (no omitted helicity factor):

```math
\begin{aligned}
A&=3M+s-2u-v, &B&=3M+s-u-2v,\\
C&=2s(M+s-u-v)-(M+s-u)(M+s-v),\\
D&=M^3-2M^2s+Ms^2-3Muv-suv+u^2v+uv^2,\\
T&=uvCAB-4(xyvA+xuzB-yzC)D,\\
H(s,u,v)&=-\frac{T}{y\sqrt{u-4M}\,z\sqrt{v-4M}\,
                         \lambda(s,u,M)\lambda(s,v,M)}.
\end{aligned}
```

On the real physical domain the repository's additive interference is

```math
I(s)=\frac{1}{2\pi(8\pi)^2s}
\int_{4M}^{b}du\int_{v_-}^{v_+}dv\,
\operatorname{Re}[f^-(u)f^+(v)]\,H(s,u,v),\qquad
\rho_{\rm SYMM}(s)=\rho_{\rm QTB}(s)+I(s).
```

The paper instead labels the opposite-signed term ``\rho_{\rm INT}``, so its convention is ``\rho_{\rm SYMM}=\rho_{\rm QTB}-\rho_{\rm INT}``. Apparent endpoint divisions in ``H`` cancel against other factors; do not evaluate the separate expressions exactly at degenerate Dalitz corners. The historical Cuba evaluator samples the interior.

## First-sheet loop, including the physical boundary

First define the once-subtracted integral

```math
\Pi_0(s)=\frac{s}{\pi}\int_a^\infty\frac{\rho(x)}{x(x-s)}dx,
\qquad
\Pi(s)=\Pi_0(s)-\operatorname{Re}\Pi_0(s_0+i0).
```

For a real argument above threshold, `selfenergy` returns the **upper rim**:

```math
\Pi(s+i0)=\frac{1}{\pi}\left[
\int_a^\infty\frac{s[\rho(x)-\rho(s)]}{x(x-s)}dx
+\rho(s)\log\frac{a}{|a-s|}\right]
-\operatorname{Re}\Pi_0(s_0+i0)+i\rho(s).
```

The integrand in brackets has removable limit ``\rho'(s)`` at ``x=s``; this is a regular quadrature, not a finite imaginary regulator. Below threshold the original real integral is regular. At threshold this implementation requires ``\rho(a)=0`` and a convergent threshold integral. The class accepts continuous, piecewise differentiable densities with growth slower than ``x`` at infinity. A cusp inside the integration interval needs appropriate quadrature care.

The package's physical-density lookup uses 100 points uniformly spaced in ``\arctan(s-a)`` up to ``10^4`` GeV², linear interpolation, zero below threshold, and a constant ``1/(8\pi)`` above the last point. This is the historical approximation, not an accuracy certificate. Tightening the quadrature tolerance cannot remove interpolation error. Use `DispersiveLineshape` with a refined density to assess that error. The artifact metadata records the default approximation.

## Amplitude, spectrum, and adjacent sheet

```math
t(s)=\frac{g^2}{m^2-s-\tfrac12g^2\Pi(s)},\qquad
\mathcal I(s)=|n t(s)|^2\rho(s)\frac{(m_\tau^2-s)^2}{s}
\left(1+\frac{2s}{m_\tau^2}\right).
```

The spectrum is zero outside ``a<s<m_\tau^2``. On the upper rim,
``\operatorname{Im}t^{-1}=-\rho/2``. The legacy functions use
``\rho_{\rm disp}=-i\Pi``; this accounts for their explicit ``i`` in the denominator.

For continuation from the upper rim into the adjacent lower-half-plane sheet,

```math
\Pi_{\rm II}(s)=\Pi_{\rm I}(s)+2i\rho_{\rm cont}(s),\qquad
 t_{\rm II}(s)=\frac{g^2}{m^2-s-\tfrac12g^2\Pi_{\rm II}(s)}.
```

Here ``\rho_{\rm cont}`` must be the holomorphic continuation with the intended contour history. The public API requires it explicitly and restricts this sheet-II formula to ``\operatorname{Im}s<0``. A square root or an interpolated complex grid alone does not certify the correct sheet. The addendum identifies remaining global transport issues. No new published pole is claimed by this modernization.

## Reusable artifacts

`src/Dispersive.jl` is a reusable generic implementation requiring only QuadGK. `scripts/export_lineshape.jl` writes a CSV with ``s,\rho,\operatorname{Re}\Pi,\operatorname{Im}\Pi,\operatorname{Re}t,\operatorname{Im}t,\mathcal I`` and a TOML metadata file. The equations on this page are the language-independent specification. `scripts/construct_JPACModel.jl` exports both QTB and SYMM examples. Neither export fits the input data: its ``m=1.3`` GeV, ``g=7`` GeV, ``n=1`` are illustrative.

A ready-to-read table and metadata are committed in [`examples/`](https://github.com/mmikhasenko/PolePositionOfA1fromTau/tree/main/examples). Regenerate them with:

```sh
A1_OUTPUT_DIR=examples julia +1.13 --project=. scripts/export_lineshape.jl
```
