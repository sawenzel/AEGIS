#!/usr/bin/env python3
"""Generate the C++ bodies of the diffcross.f integral functions.

Emits include/DiffCrossGenerated.inc, which include/DiffCross.h pulls in. The
generated file IS committed, so building needs no Python; regenerate and diff
it whenever the Fortran changes.

Division of labour, deliberately:
  * generated  -- every expression, in particular the N1..N18 polynomials
                  (N12 alone is 77 continuation lines) and the Id2/Id3 bodies.
  * hand-written in DiffCross.h -- the two places with control flow worth
                  reading (Iz0's guard and the tail of diffCross), plus the
                  scaffolding: ipow, Vec2, PhysParams, pi.

Transcription is the highest-risk step of the port: one sign or subscript slip
yields a plausible wrong cross section that survives every smoke test. So it is
mechanical and reproducible rather than careful-by-hand.

Usage:  gen_diffcross.py <diffcross.f> -o <out.inc>
"""
from __future__ import annotations

import argparse
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from f2cpp import (COMMENT_RE, convert_expr, logical_lines,  # noqa: E402
                   split_routines)

# Canonical spellings restored after identifiers are lowercased (Fortran is
# case-insensitive, C++ is not, and diffcross.f declares RX but assigns rx).
CANON = {f"i{k}{n}": f"I{k}{n}"
         for k, ns in (("z", "012"), ("d", "0123"), ("v", "012")) for n in ns}

# Supplied by the scaffolding, not declared per function.
GLOBAL_NAMES = {"pi"}

# `DOUBLE PRECISION Iz2,Id1,Id2,Id3,Iv1,Iv2` inside Diffcross declares the
# TYPE OF THE FUNCTIONS, not locals. Declaring T variables with those names
# would shadow the functions and fail to compile in a confusing way.
FUNCTION_NAMES = set(CANON)

# COMMON /PHYSPARAM/ members: become const aliases of the PhysParams argument.
PHYSPARAM = ["gamma", "beta", "m", "w1xw1", "w1xw2", "w2xw2", "wl", "wy"]

# The 18 amplitude terms, summed into NT. Their maximum magnitude is the
# denominator of the cancellation monitor: |NT| / max|Ni| measures how much
# significance the sum destroyed, which is the trigger for escalating to a
# wider type. See doc/02-findings.md section 3.
NTERMS = [f"n{i}" for i in range(1, 19)]

# Every integral function can hit Iz0's unphysical branch, so the flag is
# threaded through all of them rather than smuggled in a global.
OK_PARAM = "bool& ok"

SIG_RE = re.compile(
    r"^\s*(?:DOUBLE\s+PRECISION\s+)?(?:SUBROUTINE|FUNCTION)\s+(\w+)\s*\(([^)]*)\)",
    re.I)
# Must not match `DOUBLE PRECISION FUNCTION Iz0(x,u,v)`, or FUNCTION and the
# routine name get collected as if they were local variables.
DECL_RE = re.compile(r"^\s*DOUBLE\s+PRECISION\s+(?!FUNCTION\b)(.+)$", re.I)


def lower_idents(text: str) -> str:
    """Lowercase identifiers, leaving literals and the T template alone."""
    return re.sub(r"\b([A-Za-z]\w*)\b",
                  lambda m: m.group(1) if m.group(1) == "T" else m.group(1).lower(),
                  text)


def recanon(text: str) -> str:
    for low, hi in CANON.items():
        text = re.sub(rf"\b{low}\b", hi, text)
    # pi is a variable template in the scaffolding, so it needs its argument.
    return re.sub(r"\bpi\b", "pi<T>", text)


def parse_routine(body, name):
    """Return (args, declared) for one routine, from ITS OWN statements only.

    `body` must be the routine's statement list (signature first), not the
    whole file -- scanning the file mixes every routine's declarations
    together and every function ends up with 150 locals.
    """
    args, declared = [], {}          # declared: name -> is_array
    for _, s in body:
        m = SIG_RE.match(s)
        if m and m.group(1).lower() == name.lower():
            args = [a.strip().lower() for a in m.group(2).split(",") if a.strip()]
            continue
        m = DECL_RE.match(s)
        if m:
            # split on commas not inside parentheses
            for item in re.findall(r"(\w+\s*(?:\(\s*\d+\s*\))?)", m.group(1)):
                item = item.strip()
                if not item:
                    continue
                nm = re.match(r"(\w+)", item).group(1).lower()
                declared[nm] = "(" in item
    return args, declared


def emit_function(name, args, declared, assignments):
    arrays = {k for k, v in declared.items() if v}
    ret = name.lower()

    params = []
    for a in args:
        params.append(f"const Vec2<T>& {a}" if a in arrays else f"T {a}")
    params.append(OK_PARAM)

    locals_ = [k for k in declared
               if k not in args and k not in GLOBAL_NAMES and k != ret
               and k not in FUNCTION_NAMES]
    lines = [f"template <typename T>",
             f"T {recanon(name.lower())}({', '.join(params)}) {{"]
    vec_locals = sorted(k for k in locals_ if k in arrays)
    sca_locals = sorted(k for k in locals_ if k not in arrays)
    if vec_locals:
        lines.append("  Vec2<T> " + ", ".join(vec_locals) + ";")
    if sca_locals:
        # Zero-initialise: several routines assign only some locals on some
        # paths, and an uninitialised read is UB in C++ where Fortran merely
        # gave you whatever was on the stack.
        lines.append("  T " + ", ".join(f"{s}{{}}" for s in sca_locals) + ";")
    lines.append(f"  T {ret}{{}};")
    for lhs, rhs in assignments:
        lines.append(f"  {lhs} = {rhs};")
    lines.append(f"  return {ret};")
    lines.append("}")
    return "\n".join(lines)


def assignments_of(body):
    """(lhs, rhs) pairs, converted, lowercased and re-canonicalised."""
    skip = re.compile(
        r"^\s*(IMPLICIT|DOUBLE\s+PRECISION|INTEGER|COMMON|DIMENSION|PARAMETER|"
        r"EXTERNAL|DATA|RETURN|IF|ELSEIF|ELSE\s*IF|ELSE|ENDIF|END\s*IF|THEN|"
        r"CALL|WRITE|GO\s*TO|GOTO|CHARACTER|CONTINUE|SUBROUTINE|FUNCTION)\b",
        re.I)
    out = []
    for lineno, s in body:
        if skip.match(s):
            continue
        m = re.match(r"^\s*([A-Za-z]\w*(?:\s*\(\s*\d+\s*\))?)\s*=\s*(.+)$", s)
        if not m:
            raise ValueError(f"line {lineno}: unrecognised: {s[:80]!r}")
        # A second `=` in the RHS means two statements were glued together by
        # the line lexer -- the tab-continuation trap. Catch it here: such a
        # merge can produce code that still compiles and is simply wrong.
        if re.search(r"(?<![<>=!])=(?!=)", m.group(2)):
            raise ValueError(
                f"line {lineno}: RHS contains '=', statements were probably "
                f"merged by the continuation lexer: {s[:100]!r}")
        # LHS is deliberately NOT re-canonicalised. Fortran assigns the result
        # to the function's own name, so recanon would turn `iz1 = ...` into
        # `Iz1 = ...` -- assigning to the function rather than the result
        # variable. Declarations use the lowercase form throughout.
        lhs = convert_expr(lower_idents(m.group(1)))
        rhs = recanon(convert_expr(lower_idents(m.group(2))))
        # COMMON /SETZPARAM/ becomes the threaded flag. Fortran's `setzero=1`
        # means "this point is unphysical", i.e. ok = false.
        if lhs == "setzero":
            out.append(("ok", "true" if rhs.strip() in ("T(0e0)", "0") else "false"))
            continue
        # Calls to the integral functions must forward the flag.
        rhs = re.sub(r"\b(I[zdv]\d)\(([^()]*(?:\([^()]*\)[^()]*)*)\)",
                     lambda mm: f"{mm.group(1)}({mm.group(2)}, ok)", rhs)
        out.append((lhs, rhs))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    stmts = logical_lines(a.source)
    routines = split_routines(stmts)

    parts = [
        "// GENERATED by tools/gen_diffcross.py from ../TEPEMGEN/diffcross.f",
        "// Do not edit. Regenerate and diff instead.",
        "//",
        "// Iz0 and the diffCross scaffolding are hand-written in DiffCross.h:",
        "// they carry the only control flow worth reading.",
        "",
    ]
    for name in ("Iz1", "Iz2", "Id0", "Id1", "Id2", "Id3", "Iv0", "Iv1", "Iv2"):
        body = routines[name]
        args, declared = parse_routine(body, name)
        parts.append(emit_function(name, args, declared, assignments_of(body)))
        parts.append("")

    # ---- diffCrossTerms: the kinematics and the 18 amplitude terms ----
    dbody = assignments_of(routines["Diffcross"])
    # Drop the trailing guard (dsigma=NT; dsigma=0; badcount++). It is replaced,
    # not ported: discarding on the sign of dsigma is precisely what fails.
    cut = next(i for i, (l, _) in enumerate(dbody) if l == "dsigma")
    dbody = dbody[:cut]

    dargs, ddecl = parse_routine(routines["Diffcross"], "Diffcross")
    dargs = [x for x in dargs if x != "dsigma"]
    dlocals = [k for k in ddecl
               if k not in dargs and k not in GLOBAL_NAMES
               and k not in FUNCTION_NAMES and k not in PHYSPARAM
               and k != "dsigma"]
    darrays = {k for k, v in ddecl.items() if v}

    f = ["template <typename T>",
         "T diffCrossTerms(const PhysParams<T>& phys, T ppvt, T yp, T pmvt,",
         "                 T ym, T dphi, bool& ok, T* survival) {"]
    f.append("  // COMMON /PHYSPARAM/, set once by initDiffCross.")
    for nm in PHYSPARAM:
        f.append(f"  const T {nm} = phys.{nm};")
    vl = sorted(k for k in dlocals if k in darrays)
    sl = sorted(k for k in dlocals if k not in darrays)
    if vl:
        f.append("  Vec2<T> " + ", ".join(vl) + ";")
    if sl:
        f.append("  T " + ", ".join(f"{s}{{}}" for s in sl) + ";")
    emitted_monitor = False
    for lhs, rhs in dbody:
        f.append(f"  {lhs} = {rhs};")
        # The monitor must be taken at the moment of summation. Everything
        # after this point rescales NT (by 4/beta^2, by (2pi)^-4, by 1/4, and
        # into kbarn) without touching the Ni, so a ratio computed later would
        # fold those factors in and measure the wrong thing.
        if not emitted_monitor and lhs == "nt":
            emitted_monitor = True
            f.append("")
            f.append("  // Cancellation monitor: |NT| / max|Ni| is the fraction of the")
            f.append("  // sum that survived. Small means the result is noise. See")
            f.append("  // doc/02-findings.md section 3 -- this replaces the sign test.")
            f.append("  if (survival) {")
            f.append("    T mx{};")
            for t in NTERMS:
                f.append(f"    {{ const T a = {t} < T(0) ? -{t} : {t}; if (a > mx) mx = a; }}")
            f.append("    // Formed HERE, from the raw sum. Taking |NT| after the")
            f.append("    // normalisation factors below would divide a rescaled")
            f.append("    // numerator by an unscaled denominator and understate")
            f.append("    // survival by ~2.5e-4 -- enough to make every point look")
            f.append("    // catastrophically ill-conditioned.")
            f.append("    const T a = nt < T(0) ? -nt : nt;")
            f.append("    *survival = mx > T(0) ? a / mx : T(0);")
            f.append("  }")
            f.append("")
    if not emitted_monitor:
        raise ValueError("never found the `NT = N1+...+N18` summation")
    f.append("  return nt;")
    f.append("}")
    parts.append("\n".join(f))
    parts.append("")

    open(a.out, "w").write("\n".join(parts) + "\n")
    print(f"wrote {a.out}: 9 integral functions + diffCrossTerms "
          f"({len(dbody)} statements)")


if __name__ == "__main__":
    main()
