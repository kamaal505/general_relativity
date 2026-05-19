#!/usr/bin/env python3
"""
General Relativity Calculator
==============================
Computes, from a user-supplied metric tensor:
  - Christoffel symbols (all components, including zeros)
  - Riemann curvature tensor
  - Ricci tensor
  - Ricci scalar

Outputs are saved as .txt files in a user-named folder inside the project root.
Components may be concrete expressions (e.g. -1+2*M/r) or abstract functions
(e.g. f(r), A(r)*B(theta)) — sympy handles both symbolically.
"""

import os
import re
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
    Rewrite physicist notation  sin**2(x)  →  sin(x)**2
    so that the exponent applies to the function's return value, not the
    function object itself (which would make sympy try to call an Integer).
    Handles arbitrarily nested parentheses in the argument.
    """
    result = []
    pos = 0
    for m in _POWER_FUNC_RE.finditer(expr_str):
        func_name = m.group(1)
        exponent  = m.group(2)
        arg_start = m.end()          # character after the opening '('

        depth = 1
        i = arg_start
        while i < len(expr_str) and depth:
            if expr_str[i] == "(":
                depth += 1
            elif expr_str[i] == ")":
                depth -= 1
            i += 1
        arg_end = i - 1              # position of the matching ')'

        arg = expr_str[arg_start:arg_end]
        result.append(expr_str[pos:m.start()])
        result.append(f"{func_name}({arg})**{exponent}")
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

    # Promote unknown call-position identifiers to sympy Functions
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
    while True:
        raw = clean(input("  > ")).split()
        if len(raw) == n:
            try:
                return [sp.Symbol(name) for name in raw]
            except Exception as e:
                print(f"  Bad name: {e}. Try again.")
        else:
            print(f"  Need exactly {n} names. Got {len(raw)}. Try again.")


def get_metric(n: int, coords: list) -> sp.Matrix:
    coord_vars = {str(c): c for c in coords}
    coord_strs = [str(c) for c in coords]
    g = sp.zeros(n, n)

    print(f"\n  Enter the metric tensor components.")
    print(f"  Coordinate variables : {coord_strs}")
    print(f"  Built-in functions   : sin, cos, tan, exp, log, sqrt, pi, ...")
    print(f"  Abstract functions   : write f(r), A(r), h(r,theta), etc.")
    print(f"                         — sympy will differentiate them for you.")
    print(f"  The metric is symmetric -- enter components for i <= j only.")
    print(f"  Press Enter to leave a component as 0.\n")

    for i in range(n):
        for j in range(i, n):
            while True:
                try:
                    raw = clean(input(f"    g_{i}{j} = "))
                    if raw == "":
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
    """Γ^λ_μν = (1/2) g^λσ (∂_μ g_νσ + ∂_ν g_μσ − ∂_σ g_μν)"""
    total = n ** 3
    done = 0
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
    """R^ρ_σμν = ∂_μ Γ^ρ_νσ − ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ − Γ^ρ_νλ Γ^λ_μσ"""
    total = n ** 4
    done = 0
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
    """R_μν = R^λ_μλν  (trace over first and third indices of Riemann)"""
    Ric = [[sp.Integer(0)] * n for _ in range(n)]
    for mu in range(n):
        for nu in range(n):
            Ric[mu][nu] = simplify(sum(R[lam][mu][lam][nu] for lam in range(n)))
    return Ric


def compute_ricci_scalar(Ric: list, g_inv: sp.Matrix, n: int) -> sp.Expr:
    """R = g^μν R_μν"""
    R_scalar = sp.Integer(0)
    for mu in range(n):
        for nu in range(n):
            R_scalar += g_inv[mu, nu] * Ric[mu][nu]
    return simplify(R_scalar)


# ---------------------------------------------------------------------------
# Output formatting
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
    lines.append("  Convention:  Γ^λ_μν = (1/2) g^λσ (∂_μ g_νσ + ∂_ν g_μσ − ∂_σ g_μν)\n")
    lines.append("  All components (including zeros):\n")

    for lam in range(n):
        for mu in range(n):
            for nu in range(n):
                lines.append(f"    Γ^{lam}_{mu}{nu} = {Gamma[lam][mu][nu]}")
        lines.append("")

    return "\n".join(lines)


def format_riemann(R: list, n: int) -> str:
    lines = [section("RIEMANN CURVATURE TENSOR")]
    lines.append("  Convention:  R^ρ_σμν = ∂_μ Γ^ρ_νσ − ∂_ν Γ^ρ_μσ")
    lines.append("                        + Γ^ρ_μλ Γ^λ_νσ − Γ^ρ_νλ Γ^λ_μσ\n")
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
    lines.append("  Convention:  R_μν = R^λ_μλν  (contraction of Riemann on 1st & 3rd indices)\n")
    lines.append("  All components (including zeros):\n")

    for mu in range(n):
        for nu in range(n):
            lines.append(f"    R_{mu}{nu} = {Ric[mu][nu]}")
        lines.append("")

    return "\n".join(lines)


def format_ricci_scalar(R_scalar: sp.Expr) -> str:
    lines = [section("RICCI SCALAR")]
    lines.append("  Convention:  R = g^μν R_μν\n")
    lines.append(f"  R = {R_scalar}\n")
    return "\n".join(lines)


def save(content: str, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"  Saved: {path}")


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

    # --- output folder ---
    output_name = prompt("Name of output folder to create", "output")

    project_root = os.path.dirname(os.path.abspath(__file__))
    out_dir      = os.path.join(project_root, output_name)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n  Output will be written to: {out_dir}\n")

    # --- dimension & coordinates ---
    print(THIN)
    n      = get_int("Number of dimensions (e.g. 4 for spacetime)", 2, 9)
    coords = get_coordinates(n)

    # --- metric input ---
    print(THIN)
    g = get_metric(n, coords)

    # --- inverse metric ---
    print(f"\n  Computing inverse metric...", flush=True)
    try:
        g_inv = g.inv()
        g_inv = sp.simplify(g_inv)
    except Exception as e:
        print(f"\n  ERROR: Could not invert the metric. Is it non-degenerate?\n  {e}")
        sys.exit(1)

    # --- computations ---
    print(THIN)
    print("  Running tensor computations. This may take a few minutes for 4D metrics.")
    print()

    Gamma    = compute_christoffel(g, g_inv, coords, n)
    R        = compute_riemann(Gamma, coords, n)
    Ric      = compute_ricci_tensor(R, n)
    R_scalar = compute_ricci_scalar(Ric, g_inv, n)

    # --- format & save ---
    print(THIN)
    print("  Saving output files...\n")

    save(format_metric(g, coords, n),      os.path.join(out_dir, "metric.txt"))
    save(format_christoffel(Gamma, n),     os.path.join(out_dir, "christoffel.txt"))
    save(format_riemann(R, n),             os.path.join(out_dir, "riemann.txt"))
    save(format_ricci_tensor(Ric, n),      os.path.join(out_dir, "ricci_tensor.txt"))
    save(format_ricci_scalar(R_scalar),    os.path.join(out_dir, "ricci_scalar.txt"))

    print(f"\n{DIVIDER}")
    print("  All done! Five files written to:")
    print(f"    {out_dir}")
    print(f"    metric.txt  |  christoffel.txt  |  riemann.txt")
    print(f"    ricci_tensor.txt  |  ricci_scalar.txt")
    print(DIVIDER)
    print()


if __name__ == "__main__":
    main()
