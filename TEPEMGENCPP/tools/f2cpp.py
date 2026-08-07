#!/usr/bin/env python3
"""Mechanically translate the expression bodies of ../../TEPEMGEN/diffcross.f
into C++ expressions.

Why a script and not a careful afternoon: N12 alone is 77 continuation lines of
machine-generated polynomial, N3/N14/N16 are comparable, and a single sign or
subscript slip produces a plausible-looking wrong cross section that survives
every smoke test. Transcription is the highest-risk step of the whole port, so
it is done by a program whose output can be regenerated and diffed.

This is deliberately NOT a Fortran parser. It handles exactly the statement
forms that occur in diffcross.f, and raises on anything it does not recognise
rather than guessing. If it raises, extend it -- do not hand-edit the output.

Usage:
    f2cpp.py <diffcross.f> --routine Diffcross --emit-assignments
    f2cpp.py <diffcross.f> --list-routines
"""
from __future__ import annotations

import argparse
import re
import sys

# Names declared DIMENSION foo(2) in diffcross.f: subscripts become .x/.y
VEC2 = {
    "k1", "kd", "kx", "pmt", "ppt", "mk1", "mk1d", "mk1x", "mkd", "mkd1",
    "mkx", "mkx1", "x", "y", "z", "dyx", "dzx", "mx", "my", "mdyx",
}

# Intrinsics -> unqualified calls, so ADL picks the right overload for
# __float128 (quadmath) or a Boost multiprecision type. NEVER emit std::sqrt:
# on a Boost type that silently converts to double and discards the point.
INTRINSIC = {
    "sqrt": "sqrt", "cosh": "cosh", "sinh": "sinh", "log": "log",
    "cos": "cos", "sin": "sin", "abs": "abs",
}

CONT_RE = re.compile(r"^     [^ 0]")
COMMENT_RE = re.compile(r"^[Cc*!]")
ROUTINE_RE = re.compile(
    r"^\s*(?:DOUBLE\s+PRECISION\s+)?(SUBROUTINE|FUNCTION)\s+(\w+)", re.I)


def _classify(line: str):
    """(is_continuation, statement_text) for one fixed-form source line.

    Handles the GNU tab extension, which diffcross.f uses and which a naive
    "column 6 is non-blank" test gets exactly backwards. Line 778 is

        \\t  ...5 spaces then a TAB then...  Id3= -( ...

    A tab in the first six columns means the next character is column 7, i.e.
    the start of a NEW statement -- unless that character is a digit 1-9, which
    marks a continuation. Reading the tab as "column 6 is non-blank" silently
    glues `E11 = tepxx` onto the 40-line Id3 assignment that follows.
    """
    head = line[:6]
    tab = head.find("\t")
    if tab >= 0:
        nxt = line[tab + 1: tab + 2]
        if nxt.isdigit() and nxt != "0":
            return True, line[tab + 2:].strip()
        return False, line[tab + 1:].strip()
    if len(line) > 5 and line[5] not in (" ", "0"):
        return True, line[6:].strip()
    return False, line.strip()


def logical_lines(path: str):
    """Yield (first_lineno, joined_source) for each Fortran statement."""
    out, cur, start = [], "", 0
    for n, raw in enumerate(open(path), 1):
        line = raw.rstrip("\n")
        if not line.strip() or COMMENT_RE.match(line):
            continue
        cont, text = _classify(line)
        if cont:
            cur += " " + text
            continue
        if cur:
            out.append((start, cur))
        cur, start = text, n
    if cur:
        out.append((start, cur))
    return out


def split_routines(stmts):
    """Group statements by enclosing routine.

    The signature line is kept as the first statement of each body: callers
    need it to tell dummy arguments from locals, and scanning the whole file
    for declarations instead would mix every routine's locals together.
    """
    routines, name, body = {}, None, []
    for lineno, s in stmts:
        m = ROUTINE_RE.match(s)
        if m:
            if name:
                routines[name] = body
            name, body = m.group(2), [(lineno, s)]
            continue
        if re.match(r"^\s*END\s*$", s, re.I):
            if name:
                routines[name] = body
            name, body = None, []
            continue
        if name:
            body.append((lineno, s))
    if name:
        routines[name] = body
    return routines


def convert_powers(e: str) -> str:
    """`a**n` -> ipow<n>(a).

    Avoids std::pow entirely: it is slow, and on a user-supplied arithmetic
    type it either fails to compile or silently narrows. ipow<n> also evaluates
    the base once, where an expansion like `x*sq(sq(x))` would repeat it.
    """
    while True:
        i = e.find("**")
        if i < 0:
            return e
        m = re.match(r"\*\*\s*(\d+)", e[i:])
        if not m:
            raise ValueError(f"non-integer exponent near: {e[i:i+30]!r}")
        n = int(m.group(1))
        after = i + m.end()

        # Scan left for the base: a balanced parenthesised group, or an
        # identifier possibly followed by a call/subscript argument list.
        j = i
        while j > 0 and e[j - 1] == " ":
            j -= 1
        if j > 0 and e[j - 1] == ")":
            depth, k = 0, j
            while k > 0:
                k -= 1
                if e[k] == ")":
                    depth += 1
                elif e[k] == "(":
                    depth -= 1
                    if depth == 0:
                        break
            while k > 0 and (e[k - 1].isalnum() or e[k - 1] == "_"):
                k -= 1
            base_start = k
        else:
            k = j
            while k > 0 and (e[k - 1].isalnum() or e[k - 1] in "_."):
                k -= 1
            base_start = k
        base = e[base_start:j].strip()
        if not base:
            raise ValueError(f"cannot find base for ** near: {e[max(0,i-30):i+10]!r}")

        e = e[:base_start] + f"ipow<{n}>({base})" + e[after:]


def convert_expr(e: str) -> str:
    # Fortran double literals: 2D0 -> T(2e0), 1.97D0 -> T(1.97e0), 1D-3 -> T(1e-3)
    #
    # PRECISION TRAP: `T(3.14159...e0)` builds a *double* literal first and only
    # then converts, so anything beyond double precision is lost before T ever
    # sees it -- silently, and exactly in the quad build the port exists for.
    # diffcross.f has one such constant (PI, 25 digits) and it lives in a
    # PARAMETER, which is skipped as a declaration and supplied instead as a
    # templated constant. Refuse anything else rather than let it regress.
    def literal(m):
        mant, exp = m.group(1), int(m.group(2))
        digits = mant.replace(".", "").lstrip("0")
        if len(digits) > 17:
            raise ValueError(
                f"literal {m.group(0)!r} has {len(digits)} significant digits; "
                "emit it as a templated constant instead of inlining it.")
        # Emit an exact rational, not a decimal literal. `T(0.938)` parses
        # 0.938 as a *double* -- rounding it before T ever sees it -- whereas
        # gfortran parses 0.938D0 at the full width of T. The two then differ
        # by ~1e-17 relative, which is invisible in the double build and
        # swamps the comparison in the quad one. 0.938 is not dyadic, so no
        # number of digits fixes this; the digit-count check above does not
        # catch it. T(938)/T(1000) is exact for every T.
        frac = len(mant.split(".")[1]) if "." in mant else 0
        num = int(mant.replace(".", ""))
        den_pow = frac - exp
        if den_pow == 0:
            return f"T({num})"
        if den_pow > 0:
            return f"(T({num})/T({10 ** den_pow}))"
        return f"(T({num})*T({10 ** -den_pow}))"
    e = re.sub(r"(\d+(?:\.\d*)?)[Dd]([+-]?\d+)", literal, e)
    e = convert_powers(e)
    # vector subscripts -> members
    for v in sorted(VEC2, key=len, reverse=True):
        e = re.sub(rf"\b{v}\s*\(\s*1\s*\)", f"{v}.x", e)
        e = re.sub(rf"\b{v}\s*\(\s*2\s*\)", f"{v}.y", e)
    # intrinsics -> unqualified, ADL-friendly
    def fix_call(m):
        name = m.group(1)
        low = name.lower()
        return (INTRINSIC[low] if low in INTRINSIC else name) + "("
    e = re.sub(r"\b([A-Za-z]\w*)\s*\(", fix_call, e)
    if "**" in e:
        raise ValueError(f"leftover ** in: {e}")
    return re.sub(r"\s+", " ", e).strip()


def emit_assignments(body):
    """Yield (lhs, cpp_rhs) for plain assignments, skipping declarations."""
    skip = re.compile(
        r"^\s*(IMPLICIT|DOUBLE\s+PRECISION|INTEGER|COMMON|DIMENSION|PARAMETER|"
        r"EXTERNAL|DATA|RETURN|IF|ELSEIF|ELSE\s*IF|ELSE|ENDIF|END\s*IF|THEN|"
        r"CALL|WRITE|GO\s*TO|GOTO|CHARACTER|CONTINUE|SUBROUTINE|FUNCTION)\b",
        re.I)
    # Fortran statement function, e.g. `ARCOSH(x)=LOG(x+SQRT((x-1D0)*(x+1D0)))`.
    # Distinguishable from an array element assignment by the subscript being an
    # identifier rather than a literal integer. Reported, not silently dropped:
    # each one has to be hand-written into the C++ scaffolding.
    stmtfn = re.compile(r"^\s*([A-Za-z]\w*)\s*\(\s*[A-Za-z]\w*\s*\)\s*=")
    for lineno, s in body:
        if skip.match(s):
            continue
        m = stmtfn.match(s)
        if m:
            print(f"// STATEMENT FUNCTION at line {lineno}, write by hand: {s}",
                  file=sys.stderr)
            continue
        m = re.match(r"^\s*([A-Za-z]\w*(?:\s*\(\s*\d+\s*\))?)\s*=\s*(.+)$", s)
        if not m:
            raise ValueError(f"line {lineno}: unrecognised statement: {s[:80]!r}")
        # A second `=` in the RHS means two statements were glued together by
        # the line lexer -- the tab-continuation trap. Catch it here: such a
        # merge can produce code that still compiles and is simply wrong.
        if re.search(r"(?<![<>=!])=(?!=)", m.group(2)):
            raise ValueError(
                f"line {lineno}: RHS contains '=', statements were probably "
                f"merged by the continuation lexer: {s[:100]!r}")
        lhs, rhs = m.group(1), m.group(2)
        yield lineno, convert_expr(lhs), convert_expr(rhs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("--routine")
    ap.add_argument("--list-routines", action="store_true")
    ap.add_argument("--emit-assignments", action="store_true")
    a = ap.parse_args()

    routines = split_routines(logical_lines(a.source))
    if a.list_routines:
        for k, v in routines.items():
            print(f"{k:16s} {len(v)} statements")
        return
    if not a.routine or a.routine not in routines:
        sys.exit(f"routine not found; have: {', '.join(routines)}")
    for lineno, lhs, rhs in emit_assignments(routines[a.routine]):
        print(f"{lhs} = {rhs};")


if __name__ == "__main__":
    main()
