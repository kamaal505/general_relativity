#!/usr/bin/env python3
"""
General Relativity Calculator
==============================
Computes, from a user-supplied metric tensor:
  - Christoffel symbols (all components, including zeros)
  - Riemann curvature tensor
  - Ricci tensor
  - Ricci scalar

Outputs: .txt files (all components) + a .tex file rendered to PDF.
Components may be concrete expressions (e.g. -1+2*M/r) or abstract
functions (e.g. f(r), A(r)*B(theta)) — sympy handles both symbolically.
Greek coordinate names (theta, phi, rho, ...) render as proper LaTeX
symbols automatically.
"""

import os
import re
import shutil
import subprocess
import sys
import sympy as sp
from sympy import diff, Rational, simplify

# Force UTF-8 output so Unicode symbols display correctly on all terminals.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------

MATH_NAMESPACE = {
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
    "asin": sp.asin, "acos": sp.acos, "atan": sp.atan, "atan2": sp.atan2,
    "sinh": sp.sinh, "cosh": sp.cosh, "tanh": sp.tanh,
    "exp": sp.exp, "log": sp.log, "ln": sp.log,
    "sqrt": sp.sqrt, "Abs": sp.Abs, "sign": sp.sign,
    "pi": sp.pi, "E": sp.E,
}

# Strips ANSI escape sequences (arrow keys, colour codes, etc.) that some
# terminals inject into the input stream.
_ANSI_RE = re.compile(r"\x1b(?:\[[0-9;]*[A-Za-z]|[^\[A-Za-z])")


def clean(raw: str) -> str:
    """Remove ANSI escapes and surrounding whitespace, strip UTF-8 BOM."""
    return _ANSI_RE.sub("", raw).strip().lstrip("﻿")


_POWER_FUNC_RE = re.compile(r"\b([A-Za-z_]\w*)\s*\*\*\s*([A-Za-z0-9_]+)\s*\(")


def _rewrite_power_funcs(expr_str: str) -> str:
    """
    Rewrite physicist notation  sin**2(x)  ->  sin(x)**2
    so that the exponent applies to the function's return value, not the
    function object itself (which makes sympy try to call an Integer).
    Handles arbitrarily nested parentheses in the argument.
    """
    result = []
    pos = 0
    for m in _POWER_FUNC_RE.finditer(expr_str):
        func_name = m.group(1)
        exponent  = m.group(2)
        arg_start = m.end()

        depth = 1
        i = arg_start
        while i < len(expr_str) and depth:
            if expr_str[i] == "(":
                depth += 1
            elif expr_str[i] == ")":
                depth -= 1
            i += 1
        arg_end = i - 1

        result.append(expr_str[pos:m.start()])
        result.append(f"{func_name}({expr_str[arg_start:arg_end]})**{exponent}")
        pos = i

    result.append(expr_str[pos:])
    return "".join(result)


def sympify_smart(expr_str: str, coord_vars: dict) -> sp.Expr:
    """
    Parse expr_str as a sympy expression.
    - sin**2(x) style notation is rewritten to sin(x)**2 before parsing.
    - Identifiers in call position that are not built-ins are promoted to
      undefined sympy Functions so sympy can differentiate them symbolically.
    """
    expr_str  = _rewrite_power_funcs(expr_str)
    namespace = {**MATH_NAMESPACE, **coord_vars}

    called = set(re.findall(r"\b([A-Za-z_]\w*)\s*\(", expr_str))
    for name in called:
        if name not in namespace:
            namespace[name] = sp.Function(name)

    return sp.sympify(expr_str, locals=namespace)


def prompt(msg: str, default: str = "") -> str:
    label = f"{msg} [{default}]: " if default else f"{msg}: "
    val = clean(input(label))
    return val if val else default


def get_int(msg: str, minimum: int = 1, maximum: int = 10) -> int:
    while True:
        try:
            val = int(clean(input(f"{msg}: ")))
            if minimum <= val <= maximum:
                return val
            print(f"  Please enter a value between {minimum} and {maximum}.")
        except ValueError:
            print("  Please enter a whole number.")


def get_coordinates(n: int) -> list:
    print(f"\n  Enter the {n} coordinate names, separated by spaces.")
    print("  Example for 4D spacetime: t r theta phi")
    print("  Greek names (theta, phi, rho, sigma, ...) render as symbols in the PDF.")
    while True:
        raw = clean(input("  > ")).split()
        if len(raw) == n:
            try:
                return [sp.Symbol(name) for name in raw]
            except Exception as e:
                print(f"  Bad name: {e}. Try again.")
        else:
            print(f"  Need exactly {n} names. Got {len(raw)}. Try again.")


def get_metric(n: int, coords: list, diagonal: bool = False) -> sp.Matrix:
    coord_vars = {str(c): c for c in coords}
    coord_strs = [str(c) for c in coords]
    g = sp.zeros(n, n)

    print(f"\n  Enter the metric tensor components.")
    print(f"  Coordinate variables : {coord_strs}")
    print(f"  Built-in functions   : sin, cos, tan, exp, log, sqrt, pi, ...")
    print(f"  Abstract functions   : write f(r), A(r), h(r,theta), etc.")
    print(f"                         sympy will differentiate them symbolically.")
    print(f"  Physicist notation   : sin**2(theta) is understood correctly.")
    print(f"  Press Enter to leave a component as 0.\n")

    if diagonal:
        print(f"  Diagonal metric — entering {n} diagonal components only.\n")
        for i in range(n):
            while True:
                try:
                    raw = clean(input(f"    g_{i}{i} = "))
                    if not raw:
                        raw = "0"
                    g[i, i] = sympify_smart(raw, coord_vars)
                    break
                except Exception as e:
                    print(f"    Could not parse '{raw}': {e}. Try again.")
    else:
        print(f"  Symmetric metric — enter components for i <= j only.\n")
        for i in range(n):
            for j in range(i, n):
                while True:
                    try:
                        raw = clean(input(f"    g_{i}{j} = "))
                        if not raw:
                            raw = "0"
                        expr = sympify_smart(raw, coord_vars)
                        g[i, j] = expr
                        g[j, i] = expr
                        break
                    except Exception as e:
                        print(f"    Could not parse '{raw}': {e}. Try again.")
    return g


# ---------------------------------------------------------------------------
# Tensor computations
# ---------------------------------------------------------------------------

def compute_christoffel(g: sp.Matrix, g_inv: sp.Matrix, coords: list, n: int) -> list:
    """Gamma^lam_mu_nu = (1/2) g^lam_sig (d_mu g_nu_sig + d_nu g_mu_sig - d_sig g_mu_nu)"""
    total = n ** 3
    done  = 0
    Gamma = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]

    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                s = sp.Integer(0)
                for sig in range(n):
                    s += g_inv[lam, sig] * (
                        diff(g[nu, sig], coords[mu])
                        + diff(g[mu, sig], coords[nu])
                        - diff(g[mu, nu], coords[sig])
                    )
                Gamma[lam][mu][nu] = simplify(Rational(1, 2) * s)
                done += 1
                print(f"    Christoffel: {done}/{total} components done", end="\r", flush=True)

    print()
    return Gamma


def compute_riemann(Gamma: list, coords: list, n: int) -> list:
    """R^rho_sig_mu_nu = d_mu Gamma^rho_nu_sig - d_nu Gamma^rho_mu_sig + ..."""
    total = n ** 4
    done  = 0
    R = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]

    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    t1 = diff(Gamma[rho][nu][sig], coords[mu])
                    t2 = diff(Gamma[rho][mu][sig], coords[nu])
                    t3 = sum(Gamma[rho][mu][lam] * Gamma[lam][nu][sig] for lam in range(n))
                    t4 = sum(Gamma[rho][nu][lam] * Gamma[lam][mu][sig] for lam in range(n))
                    R[rho][sig][mu][nu] = simplify(t1 - t2 + t3 - t4)
                    done += 1
                    print(f"    Riemann: {done}/{total} components done", end="\r", flush=True)

    print()
    return R


def compute_ricci_tensor(R: list, n: int) -> list:
    """R_mu_nu = R^lam_mu_lam_nu"""
    Ric = [[sp.Integer(0)] * n for _ in range(n)]
    for mu in range(n):
        for nu in range(n):
            Ric[mu][nu] = simplify(sum(R[lam][mu][lam][nu] for lam in range(n)))
    return Ric


def compute_ricci_scalar(Ric: list, g_inv: sp.Matrix, n: int) -> sp.Expr:
    """R = g^mu_nu R_mu_nu"""
    R_scalar = sp.Integer(0)
    for mu in range(n):
        for nu in range(n):
            R_scalar += g_inv[mu, nu] * Ric[mu][nu]
    return simplify(R_scalar)


# ---------------------------------------------------------------------------
# Plain-text output
# ---------------------------------------------------------------------------

DIVIDER = "=" * 64
THIN    = "-" * 64


def section(title: str) -> str:
    return f"\n{DIVIDER}\n  {title}\n{DIVIDER}\n"


def format_metric(g: sp.Matrix, coords: list, n: int) -> str:
    lines = [section("INPUT METRIC")]
    lines.append(f"  Coordinates : {[str(c) for c in coords]}\n")
    lines.append("  g_ij components:\n")
    for i in range(n):
        for j in range(n):
            lines.append(f"    g_{i}{j} = {g[i, j]}")
    return "\n".join(lines)


def format_christoffel(Gamma: list, n: int) -> str:
    lines = [section("CHRISTOFFEL SYMBOLS")]
    lines.append("  Convention:  Gamma^lam_mu_nu = (1/2) g^lam_sig (d_mu g_nu_sig + d_nu g_mu_sig - d_sig g_mu_nu)\n")
    lines.append("  All components (including zeros):\n")
    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                lines.append(f"    Gamma^{lam}_{mu}{nu} = {Gamma[lam][mu][nu]}")
        lines.append("")
    return "\n".join(lines)


def format_riemann(R: list, n: int) -> str:
    lines = [section("RIEMANN CURVATURE TENSOR")]
    lines.append("  Convention:  R^rho_sig_mu_nu = d_mu Gamma^rho_nu_sig - d_nu Gamma^rho_mu_sig + ...\n")
    lines.append("  All components (including zeros):\n")
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    lines.append(f"    R^{rho}_{sig}{mu}{nu} = {R[rho][sig][mu][nu]}")
            lines.append("")
        lines.append("")
    return "\n".join(lines)


def format_ricci_tensor(Ric: list, n: int) -> str:
    lines = [section("RICCI TENSOR")]
    lines.append("  Convention:  R_mu_nu = R^lam_mu_lam_nu\n")
    lines.append("  All components (including zeros):\n")
    for mu in range(n):
        for nu in range(n):
            lines.append(f"    R_{mu}{nu} = {Ric[mu][nu]}")
        lines.append("")
    return "\n".join(lines)


def format_ricci_scalar(R_scalar: sp.Expr) -> str:
    lines = [section("RICCI SCALAR")]
    lines.append("  Convention:  R = g^mu_nu R_mu_nu\n")
    lines.append(f"  R = {R_scalar}\n")
    return "\n".join(lines)


def save(content: str, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# LaTeX / PDF output
# ---------------------------------------------------------------------------

def _lx(expr) -> str:
    """Shorthand: sympy expression -> LaTeX string."""
    return sp.latex(expr)


def _cidx(coord) -> str:
    """LaTeX string for a coordinate used as a tensor index.
    sp.latex automatically converts theta -> \\theta, phi -> \\phi, etc."""
    return _lx(coord)


def _dmath(lhs: str, rhs_expr) -> list:
    """Return lines for a numbered equation* block."""
    return [r"\begin{equation*}", "  " + lhs + " = " + _lx(rhs_expr), r"\end{equation*}"]


def generate_latex(
    g: sp.Matrix,
    Gamma: list,
    R: list,
    Ric: list,
    R_scalar: sp.Expr,
    coords: list,
    n: int,
) -> str:
    """Build a complete LaTeX document string for all computed quantities."""

    def ci(i: int) -> str:
        return "{" + _cidx(coords[i]) + "}"

    L = []

    # ------ preamble ------
    L += [
        r"\documentclass[12pt]{article}",
        r"\usepackage{amsmath, amssymb}",
        r"\usepackage[margin=1in]{geometry}",
        r"\setlength{\parindent}{0pt}",
        r"\setlength{\parskip}{6pt}",
        r"\title{\textbf{General Relativity: Tensor Calculations}}",
        r"\date{}",
        r"\begin{document}",
        r"\maketitle",
    ]

    coord_str = r",\ ".join("$" + _lx(c) + "$" for c in coords)
    L.append(r"\textbf{Coordinates:} " + coord_str + r"\par")

    # ------ metric ------
    L += [
        r"\section{Metric Tensor}",
        r"\begin{equation*}",
        r"  g_{\mu\nu} = " + _lx(g),
        r"\end{equation*}",
        r"\textbf{Non-zero components:}",
    ]

    metric_entries = [
        "  g_{" + ci(i) + ci(j) + "} &= " + _lx(g[i, j])
        for i in range(n) for j in range(i, n) if g[i, j] != 0
    ]
    if metric_entries:
        L.append(r"\begin{align*}")
        L.append(r" \\".join(metric_entries))
        L.append(r"\end{align*}")
    else:
        L.append(r"\textit{All metric components are zero.}")

    # ------ Christoffel ------
    L += [
        r"\section{Christoffel Symbols}",
        r"\textbf{Convention:}",
        r"\begin{equation*}",
        (r"  \Gamma^{\lambda}{}_{\mu\nu} = \frac{1}{2}\,g^{\lambda\sigma}"
         r"\!\left(\partial_\mu g_{\nu\sigma} + \partial_\nu g_{\mu\sigma}"
         r" - \partial_\sigma g_{\mu\nu}\right)"),
        r"\end{equation*}",
        r"\textbf{Non-zero components:}",
    ]

    any_nz = False
    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                val = Gamma[lam][mu][nu]
                if val != 0:
                    any_nz = True
                    lhs = r"\Gamma^{" + ci(lam) + r"}{}_{" + ci(mu) + ci(nu) + "}"
                    L += _dmath(lhs, val)

    L.append(
        r"\textit{All remaining components are zero.}"
        if any_nz else
        r"\textit{All Christoffel symbols vanish (flat spacetime).}"
    )

    # ------ Riemann ------
    L += [
        r"\section{Riemann Curvature Tensor}",
        r"\textbf{Convention:}",
        r"\begin{equation*}",
        (r"  R^{\rho}{}_{\sigma\mu\nu} ="
         r"  \partial_\mu\Gamma^{\rho}{}_{\nu\sigma}"
         r" - \partial_\nu\Gamma^{\rho}{}_{\mu\sigma}"
         r" + \Gamma^{\rho}{}_{\mu\lambda}\Gamma^{\lambda}{}_{\nu\sigma}"
         r" - \Gamma^{\rho}{}_{\nu\lambda}\Gamma^{\lambda}{}_{\mu\sigma}"),
        r"\end{equation*}",
        r"\textbf{Non-zero components:}",
    ]

    any_nz = False
    for rho in range(n):
        for sig in range(n):
            for mu in range(n):
                for nu in range(n):
                    val = R[rho][sig][mu][nu]
                    if val != 0:
                        any_nz = True
                        lhs = "R^{" + ci(rho) + "}{}_{" + ci(sig) + ci(mu) + ci(nu) + "}"
                        L += _dmath(lhs, val)

    L.append(
        r"\textit{All remaining components are zero.}"
        if any_nz else
        r"\textit{All Riemann tensor components vanish (flat spacetime).}"
    )

    # ------ Ricci tensor ------
    L += [
        r"\section{Ricci Tensor}",
        r"\textbf{Convention:} $R_{\mu\nu} = R^{\lambda}{}_{\mu\lambda\nu}$",
        r"\textbf{Non-zero components:}",
    ]

    any_nz = False
    for mu in range(n):
        for nu in range(n):
            val = Ric[mu][nu]
            if val != 0:
                any_nz = True
                lhs = "R_{" + ci(mu) + ci(nu) + "}"
                L += _dmath(lhs, val)

    L.append(
        r"\textit{All remaining components are zero.}"
        if any_nz else
        r"\textit{All Ricci tensor components vanish.}"
    )

    # ------ Ricci scalar ------
    L += [
        r"\section{Ricci Scalar}",
        r"\textbf{Convention:} $R = g^{\mu\nu} R_{\mu\nu}$",
    ]
    L += _dmath("R", R_scalar)

    L.append(r"\end{document}")
    return "\n".join(L)


def compile_pdf(tex_path: str, out_dir: str) -> bool:
    """Try to compile tex_path with pdflatex. Clean up auxiliary files."""
    if not shutil.which("pdflatex"):
        print("  pdflatex not found on PATH — .tex file saved for manual compilation.")
        print("  Tip: paste the .tex file into overleaf.com for an instant PDF.")
        return False

    print("  Compiling PDF with pdflatex ...", flush=True)
    proc = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-output-directory", out_dir, tex_path],
        capture_output=True,
        text=True,
        cwd=out_dir,
    )

    stem = os.path.splitext(os.path.basename(tex_path))[0]
    for ext in (".aux", ".log", ".out"):
        aux = os.path.join(out_dir, stem + ext)
        if os.path.exists(aux):
            os.remove(aux)

    if proc.returncode == 0:
        print(f"  PDF ready: {os.path.join(out_dir, stem + '.pdf')}")
        return True

    print("  pdflatex reported errors — the .tex file is saved for inspection.")
    for line in proc.stdout.strip().splitlines()[-10:]:
        if line.strip():
            print(f"    {line}")
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print()
    print(DIVIDER)
    print("  General Relativity Calculator")
    print("  Christoffel  |  Riemann  |  Ricci Tensor  |  Ricci Scalar")
    print(DIVIDER)
    print()

    output_name  = prompt("Name of output folder to create", "output")
    project_root = os.path.dirname(os.path.abspath(__file__))
    out_dir      = os.path.join(project_root, output_name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n  Output will be written to: {out_dir}\n")

    print(THIN)
    n      = get_int("Number of dimensions (e.g. 4 for spacetime)", 2, 9)
    coords = get_coordinates(n)

    print(THIN)
    diagonal = prompt("Is this metric diagonal? (y/n)", "n").lower() in ("y", "yes")
    g = get_metric(n, coords, diagonal=diagonal)

    print(f"\n  Computing inverse metric...", flush=True)
    try:
        g_inv = simplify(g.inv())
    except Exception as e:
        print(f"\n  ERROR: Could not invert the metric. Is it non-degenerate?\n  {e}")
        sys.exit(1)

    print(THIN)
    print("  Running tensor computations. This may take a few minutes for 4D metrics.")
    print()

    Gamma    = compute_christoffel(g, g_inv, coords, n)
    R        = compute_riemann(Gamma, coords, n)
    Ric      = compute_ricci_tensor(R, n)
    R_scalar = compute_ricci_scalar(Ric, g_inv, n)

    print(THIN)
    print("  Saving output files...\n")

    save(format_metric(g, coords, n),   os.path.join(out_dir, "metric.txt"))
    save(format_christoffel(Gamma, n),  os.path.join(out_dir, "christoffel.txt"))
    save(format_riemann(R, n),          os.path.join(out_dir, "riemann.txt"))
    save(format_ricci_tensor(Ric, n),   os.path.join(out_dir, "ricci_tensor.txt"))
    save(format_ricci_scalar(R_scalar), os.path.join(out_dir, "ricci_scalar.txt"))

    tex_path = os.path.join(out_dir, "results.tex")
    save(generate_latex(g, Gamma, R, Ric, R_scalar, coords, n), tex_path)

    print()
    compile_pdf(tex_path, out_dir)

    print(f"\n{DIVIDER}")
    print("  All done! Files written to:")
    print(f"    {out_dir}")
    print(f"    .txt: metric | christoffel | riemann | ricci_tensor | ricci_scalar")
    print(f"    .tex: results.tex  (+ results.pdf if pdflatex was found)")
    print(DIVIDER)
    print()


if __name__ == "__main__":
    main()
