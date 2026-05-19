# General Relativity Calculator

A symbolic GR calculator that takes a metric tensor as input and computes:

- Christoffel symbols $\Gamma^\lambda{}_{\mu\nu}$
- Riemann curvature tensor $R^\rho{}_{\sigma\mu\nu}$
- Ricci tensor $R_{\mu\nu}$
- Ricci scalar $R$

Results are saved as plain `.txt` files (all components, including zeros) and as a compiled **PDF** with proper tensor notation and rendered mathematics.

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
