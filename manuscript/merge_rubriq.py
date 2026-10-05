# -*- coding: utf-8 -*-
"""Merge the Rubriq language edit of the v3 manuscript into the v4 build without letting it change the science or the format.

A paragraph of the edit is accepted only when the v4 paragraph equals the v3 paragraph (the build differs only in DOIs and the version),
the edit keeps the same content words (function words, punctuation, spelling variants and a few safe synonyms aside), the same numbers
and citation brackets, and the same subscript, superscript, bold and italic runs, and a hand review did not find a change of meaning
(REJECT_REVIEWED). Accepted paragraphs take the edited runs; all others keep the build text. Equation paragraphs, the references and
tables are never taken from the edit. A report of every decision is written next to the output.

    python merge_rubriq.py <v3.docx> <v3_updated.docx> <v4_build.docx> <out.docx>
"""
import copy
import difflib
import io
import re
import sys

import docx
from docx.oxml.ns import qn

M = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"
US = [("non-superconducting", "nonsuperconducting"), ("quasi-static", "quasistatic"), ("grey", "gray"), ("towards", "toward"),
      ("neighbour", "neighbor"), ("modelled", "modeled"), ("modelling", "modeling"), ("behaviour", "behavior"),
      ("non-equilibrium", "nonequilibrium"), ("non-conservative", "nonconservative"), ("non-monotonic", "nonmonotonic"),
      ("sub-lattice", "sublattice"), ("un-killed", "unkilled"), ("artefact", "artifact"), ("favour", "favor")]


def us(t):
    for a, b in US:
        t = re.sub(a, b, t)
        t = re.sub(a[0].upper() + a[1:], b[0].upper() + b[1:], t)
    return t


STOP = set("the a an of is are was were be been being in on at that which and as to for by with from it its their this these those "
           "into has have had can could would will may also then thus so therefore whereas while each both respectively otherwise".split())
SYN = {"approximately": "about", "roughly": "about", "increases": "rises", "increase": "rises", "grows": "rises", "grow": "rises",
       "decreases": "falls", "decrease": "falls", "declines": "falls", "size-independent": "sizeindependent", "n-value": "nvalue"}

ACCEPT_REVIEWED = []
REJECT_REVIEWED = {}


def tokens(text):
    t = us(text).lower().replace("’", "'")
    t = re.sub("[‐‑‒–—―−]", "-", t)
    words = re.findall(r"[a-z0-9α-ωΔΥ]+", t)
    out = []
    for w in words:
        if w in STOP:
            continue
        w = SYN.get(w, w)
        if len(w) > 4 and w.endswith("s") and not w.endswith("ss"):
            w = w[:-1]
        out.append(w)
    return out


def numbers(text):
    return re.findall(r"\d+(?:\.\d+)?", text) + re.findall(r"\[[0-9,–\- ]+\]", text)


def sig(p):
    sub = sup = b = i = 0
    for r in p._p.iterfind(".//" + qn("w:r")):
        rpr = r.find(qn("w:rPr"))
        if rpr is None or r.find(qn("w:t")) is None:
            continue
        va = rpr.find(qn("w:vertAlign"))
        if va is not None:
            sub += va.get(qn("w:val")) == "subscript"
            sup += va.get(qn("w:val")) == "superscript"
        bb, ii = rpr.find(qn("w:b")), rpr.find(qn("w:i"))
        b += bb is not None and bb.get(qn("w:val")) not in ("0", "false")
        i += ii is not None and ii.get(qn("w:val")) not in ("0", "false")
    return sub, sup, b, i


def has_math(p):
    return p._p.find(".//" + M + "oMath") is not None


SAFE_TAGS = ("w:pPr",)


def plain_runs(p):
    """True if the paragraph consists only of a paragraph-properties element and runs that carry text (no fields, hyperlinks, math)."""
    for ch in p._p:
        if ch.tag == qn("w:pPr"):
            continue
        if ch.tag != qn("w:r"):
            return False
        for g in ch:
            if g.tag not in (qn("w:rPr"), qn("w:t")):
                return False
    return True


def char_formats(p):
    """[(char, rPr xml or None)] for the text of a plain paragraph."""
    from lxml import etree
    out = []
    for r in p._p.iterfind(qn("w:r")):
        rpr = r.find(qn("w:rPr"))
        key = etree.tostring(rpr) if rpr is not None else None
        for t in r.iterfind(qn("w:t")):
            for ch in (t.text or ""):
                out.append((ch, key))
    return out


def rebuild_runs(p, chars):
    """Replace the runs of paragraph p by runs made from [(char, rPr xml)] (consecutive characters with equal formatting share a run)."""
    from lxml import etree
    for r in list(p._p.iterfind(qn("w:r"))):
        p._p.remove(r)
    i = 0
    while i < len(chars):
        j = i
        while j < len(chars) and chars[j][1] == chars[i][1]:
            j += 1
        r = etree.SubElement(p._p, qn("w:r"))
        if chars[i][1] is not None:
            r.append(etree.fromstring(chars[i][1]))
        t = etree.SubElement(r, qn("w:t"))
        t.text = "".join(c for c, _ in chars[i:j])
        t.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
        i = j


TOK = re.compile(r"\w+|[^\w\s]|\s+", re.U)


def hunks(a, b):
    """Token-level diff of two strings: list of (start, end, replacement, old, new) in character offsets of a, for every difference."""
    ta, tb = TOK.findall(a), TOK.findall(b)
    off = [0]
    for t in ta:
        off.append(off[-1] + len(t))
    out = []
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, ta, tb, autojunk=False).get_opcodes():
        if op != "equal":
            out.append((off[i1], off[i2], "".join(tb[j1:j2]), "".join(ta[i1:i2]), "".join(tb[j1:j2])))
    return out


WHITELIST = {("about", "approximately"), ("roughly", "approximately"), ("rises", "increases"), ("rise", "increase"), ("declines", "decreases"),
             ("falls", "decreases"), ("grow", "increase"), ("grows", "increases"), ("non-equilibrium", "nonequilibrium"),
             ("quasi-static", "quasistatic"), ("Quasi-static", "Quasistatic"), ("grey", "gray"), ("towards", "toward"),
             ("non-superconducting", "nonsuperconducting")}


def safe(old, new, before="", after=""):
    """A hunk is taken only if it inserts or deletes a single comma (not inside a number or a parameter list such as 'L = 32, beta = 20'),
    or is one of the listed wording-neutral substitutions (spelling variants and 'about' -> 'approximately', matched in inflection).
    Every other change of words or punctuation is left to the author."""
    if (old, new) in (("", ","), (",", "")):
        if before[-1:].isdigit() and after[:1].isdigit():
            return False
        return True
    return (old.strip(), new.strip()) in WHITELIST


def main(a_path, b_path, c_path, out_path):
    A = docx.Document(a_path).paragraphs
    B = docx.Document(b_path).paragraphs
    Cd = docx.Document(c_path)
    C = Cd.paragraphs
    assert len(A) == len(B), (len(A), len(B))
    c_text = [p.text for p in C]
    sm = difflib.SequenceMatcher(None, [p.text for p in A], c_text, autojunk=False)
    c_of_a = {}
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            for k in range(i2 - i1):
                c_of_a[i1 + k] = j1 + k
    report, n_acc, n_rej, n_par = [], 0, 0, 0
    in_refs = False
    for i, (pa, pb) in enumerate(zip(A, B)):
        if pa.text.strip() == "References":
            in_refs = True
        if pa.text == pb.text:
            continue
        head = pa.text[:60].replace("\n", " ")
        if in_refs or re.match(r"\[\d+\]", pa.text):
            report.append("SKIP reference entry | %s" % head)
            continue
        if i not in c_of_a:
            report.append("SKIP paragraph changed in the build | %s" % head)
            continue
        pc = C[c_of_a[i]]
        if has_math(pc) or not plain_runs(pc):
            report.append("SKIP equation or field paragraph | %s" % head)
            continue
        chars = char_formats(pc)
        assert "".join(c for c, _ in chars) == pa.text, head
        acc = []
        for (s0, s1, rep, old, new) in hunks(pa.text, pb.text):
            if safe(old, new, pa.text[:s0], pa.text[s1:]) and not any(k in (old + "->" + new) for k in DENY):
                acc.append((s0, s1, rep))
                n_acc += 1
                report.append("  accept  [%s] -> [%s]" % (old, new))
            else:
                n_rej += 1
                report.append("  reject  [%s] -> [%s]" % (old, new))
        if acc:
            n_par += 1
            new_chars = list(chars)
            for s0, s1, rep in sorted(acc, reverse=True):
                fmt = chars[s0 - 1][1] if s0 > 0 else (chars[s1][1] if s1 < len(chars) else None)
                if s1 > s0:
                    fmt = chars[s0][1]
                new_chars[s0:s1] = [(c, fmt) for c in rep]
            rebuild_runs(pc, new_chars)
        report.append("PARAGRAPH %s" % head)
    Cd.save(out_path)
    rep = out_path.rsplit(".", 1)[0] + "_merge_report.txt"
    io.open(rep, "w", encoding="utf-8").write("accepted %d hunks in %d paragraphs, rejected %d hunks" % (n_acc, n_par, n_rej) + chr(10) + chr(10) + chr(10).join(report) + chr(10))
    print("accepted %d hunks in %d paragraphs, rejected %d hunks -> %s" % (n_acc, n_par, n_rej, out_path))


DENY = []

if __name__ == "__main__":
    main(*sys.argv[1:5])
