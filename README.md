# General Relativity Calculator

A symbolic GR calculator that takes a metric tensor as input and computes:

- Christoffel symbols $\Gamma^\lambda{}_{\mu\nu}$
- Riemann curvature tensor $R^\rho{}_{\sigma\mu\nu}$
- Ricci tensor $R_{\mu\nu}$
- Ricci scalar $R$
- Einstein tensor $G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}R\,g_{\mu\nu}$

Results are saved as plain `.txt` files (all components, including zeros) and as a compiled **PDF** with proper tensor notation and rendered mathematics.

Built with black hole research in mind: charged (Reissner–Nordström, Kerr–Newman), rotating (Kerr), and other non-diagonal, multi-parameter metrics are all supported — see [Choosing constant names](#choosing-constant-names) below for the one gotcha to know about.

---

## Requirements

**Python** — install dependencies with:

```bash
cd general_relativity
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt   # Windows
# or
.venv/bin/pip install -r requirements.txt       # Mac/Linux
```

**PDF compilation** (optional) — install [MiKTeX](https://miktex.org/) or [TeX Live](https://tug.org/texlive/) so that `pdflatex` is on your PATH. If it isn't found, the `.tex` file is saved and can be compiled manually or pasted into [Overleaf](https://overleaf.com).

---

## Usage

```bash
.venv\Scripts\python gr_calculator.py    # Windows
.venv/bin/python gr_calculator.py        # Mac/Linux
```

The script will ask for:

| Prompt | Example |
|--------|---------|
| Output folder name | `schwarzschild` |
| Number of dimensions | `4` |
| Coordinate names | `t r theta phi` |
| Is the metric diagonal? | `y` |
| Metric components | see below |

---

## Entering metric components

### Concrete expressions

```
g_00 = -1 + 2*M/r
g_11 = 1/(1 - 2*M/r)
```

`M` is treated as a symbolic constant automatically.

### Abstract functions

```
g_00 = f(r)
g_11 = h(r, theta)
```

Sympy differentiates abstract functions symbolically, so derivatives like $f'(r)$ and $f''(r)$ appear correctly in the output.

### Physicist power notation

```
g_33 = r**2 * sin**2(theta)
```

`sin**2(theta)` is interpreted as $\sin^2(\theta)$, not as applying `**2` to the function object.

### Supported functions

`sin cos tan sec csc cot`, their inverses (`asin`/`arcsin`, `asec`, ...), the hyperbolic set (`sinh cosh tanh sech csch coth` and inverses), plus `exp log ln sqrt Abs sign pi E`.

---

## Choosing constant names

Any name you type that isn't a coordinate or a function above is treated as a free symbolic constant — `M`, `Q`, `a`, `J`, `Lambda`, whatever your metric needs. This covers the common black hole parameters directly:

| Parameter | Typical name |
|-----------|--------------|
| Mass | `M` |
| Charge (Reissner–Nordström, Kerr–Newman) | `Q` |
| Spin (Kerr, Kerr–Newman) | `a` |
| Cosmological constant | `Lambda` |
| NUT charge, other | anything you like |

**One gotcha:** sympy itself reserves a handful of short names as built-in objects (`Q`, `I`, `N`, `O`, `S`, `C`, ...). Earlier versions of this calculator sympified those names straight to sympy's internal objects instead of a free symbol — so a charged-black-hole metric using `Q**2/r**2` would silently fail. This is now fixed: every bare name you type (that isn't one of the functions above or a coordinate) is explicitly pinned to a fresh `Symbol` before parsing, so `Q`, `I`, `N`, `S`, etc. all behave as ordinary constants. The one exception is `E`, which is kept as Euler's number (consistent with `pi`) — use a different name (e.g. `En`) if you need a constant literally called `E`.

---

## Supported symbols

Greek coordinate names are rendered as proper LaTeX symbols in the PDF automatically — no manual escaping needed.

| You type | Renders as |
|----------|-----------|
| `theta`  | $\theta$  |
| `phi`    | $\phi$    |
| `rho`    | $\rho$    |
| `sigma`  | $\sigma$  |
| `alpha`  | $\alpha$  |
| `beta`   | $\beta$   |
| `gamma`  | $\gamma$  |
| `tau`    | $\tau$    |
| `omega`  | $\omega$  |
| `mu`     | $\mu$     |
| `nu`     | $\nu$     |
| `lambda` | $\lambda$ |
| `t`, `r` | $t$, $r$  |

All other standard Greek names work too — sympy handles the full alphabet.

---

## Output files

All files are written to `general_relativity/<output_folder>/`:

| File | Contents |
|------|----------|
| `metric.txt` | Input metric components |
| `christoffel.txt` | All $\Gamma^\lambda{}_{\mu\nu}$ (including zeros) |
| `riemann.txt` | All $R^\rho{}_{\sigma\mu\nu}$ (including zeros) |
| `ricci_tensor.txt` | All $R_{\mu\nu}$ (including zeros) |
| `ricci_scalar.txt` | Ricci scalar $R$ |
| `einstein_tensor.txt` | All $G_{\mu\nu} = R_{\mu\nu} - \frac12 R g_{\mu\nu}$ (including zeros) |
| `results.tex` | Full LaTeX source |
| `results.pdf` | Compiled PDF (if pdflatex is available) |

The PDF lists only **non-zero** components with coordinate-name indices (e.g. $\Gamma^t{}_{tr}$, $R_{\theta\theta}$). The `.txt` files contain every component for completeness.

---

## Example — FLRW metric

Coordinates: `t r theta phi`, diagonal: `y`

```
g_00 = -1
g_11 = a(t)**2
g_22 = a(t)**2 * r**2
g_33 = a(t)**2 * r**2 * sin**2(theta)
```

The Ricci scalar comes out as:

$$R = \frac{6\left(a(t)\,\ddot{a}(t) + \dot{a}(t)^2\right)}{a(t)^2}$$

which is the standard FLRW result.

---

## Example — Reissner–Nordström (charged black hole)

Coordinates: `t r theta phi`, diagonal: `y`

```
g_00 = -(1 - 2*M/r + Q**2/r**2)
g_11 = 1/(1 - 2*M/r + Q**2/r**2)
g_22 = r**2
g_33 = r**2 * sin**2(theta)
```

The Ricci scalar vanishes ($R=0$, as expected — the Maxwell stress-energy tensor is traceless), while the Einstein tensor comes out charge-sourced and non-zero, e.g.

$$G_{tt} = \frac{Q^2\left(r^2 - 2Mr + Q^2\right)}{r^6}.$$

## Example — Kerr (rotating black hole)

Coordinates: `t r theta phi`, diagonal: `n` (Boyer–Lindquist form has an off-diagonal $g_{t\phi}$ term)

```
g_00 = -(1 - 2*M*r/(r**2 + a**2*cos**2(theta)))
g_03 = -2*M*r*a*sin**2(theta)/(r**2 + a**2*cos**2(theta))
g_11 = (r**2 + a**2*cos**2(theta))/(r**2 - 2*M*r + a**2)
g_22 = r**2 + a**2*cos**2(theta)
g_33 = (r**2 + a**2 + 2*M*r*a**2*sin**2(theta)/(r**2 + a**2*cos**2(theta))) * sin**2(theta)
```

This is vacuum ($R_{\mu\nu}=0$, $G_{\mu\nu}=0$), but the Riemann tensor is dense in $(r,\theta)$ — expect the computation to take noticeably longer than the diagonal examples above, since `simplify()` is called on every one of the 256 Riemann components with no symmetry shortcuts.
