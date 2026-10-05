# -*- coding: utf-8 -*-
"""Build the manuscript (docx) from the study results.

    python build_manuscript.py      ->  out/Damage_Repair_Phase_Stiffness_XY_v1.docx

Every number in Sections 2-4 is read from ../results/results.json (analyze*.py) or from the raw runs; none is typed here.
"""
import json
import math
import os
import re
import sys

import numpy as np
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "code"))
import docx_helpers as H  # noqa: E402
import meanfield as MF  # noqa: E402
import refs as RF  # noqa: E402
from docx_helpers import (DELIM, FRAC, SUB, SUP, T, V, Vs, equation, new_document)  # noqa: E402

RESDIR = os.path.join(HERE, "..", "results")
FIGDIR = os.path.join(HERE, "..", "figures")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
R = json.load(open(os.path.join(RESDIR, "results.json")))
ZEN = json.load(open(os.path.join("C:" + os.sep, "YouTube", "_dxy_zenodo_state.json")))
SW_DOI, PP_DOI = ZEN.get("software_1.1.0", ZEN["software"])["doi"], ZEN.get("publication_v2", ZEN["publication"])["doi"]
REPO = "https://github.com/sandlerleon/damage-repair-xy-stiffness"
DB = RF.load()

doc = new_document(size=11, line=1.5)
CITE = []


def cites(text):
    """[[a,b]] -> [n, m] numbered by first appearance."""
    def rep(m):
        keys = m.group(1).split(",")
        for k in keys:
            if k not in CITE:
                CITE.append(k)
        nums = sorted(CITE.index(k) + 1 for k in keys)
        out, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            out.append("%d–%d" % (nums[i], nums[j]) if j - i >= 2 else ", ".join(str(n) for n in nums[i:j + 1]))
            i = j + 1
        return "[" + ", ".join(out) + "]"
    return re.sub(r"\[\[([A-Za-z0-9_,]+)\]\]", rep, text)


LABFILE = os.path.join(OUT, "labels.json")
TLAB = json.load(open(LABFILE)) if os.path.exists(LABFILE) else {}
TLAB_NEW = {}


def tref(text):
    """@T:label@ -> 'Table n' (numbers from the previous pass; run the build twice)."""
    return re.sub(r"@T:(\w+)@", lambda m: "Table %s" % TLAB.get(m.group(1), "?"), text).replace("%%", "%")


def P(text, **kw):
    kw.setdefault("align", "justify")
    return H.para(doc, cites(tref(text)), **kw)


def HD(text, level=1):
    return H.heading(doc, text, level)


def EQ(nodes, n):
    equation(doc, nodes, n)


FIGN = [0]
TABN = [0]


def FIG(path, caption):
    FIGN[0] += 1
    H.figure(doc, os.path.join(FIGDIR, path), width_in=6.5, cap="**Figure %d.** %s" % (FIGN[0], cites(tref(caption))))


def TAB(rows, caption, widths=None, size=8.5, label=None):
    TABN[0] += 1
    if label:
        TLAB_NEW[label] = TABN[0]
    c = H.caption(doc, "**Table %d.** %s" % (TABN[0], cites(tref(caption))), keep_next=True)
    H.table(doc, rows, widths=widths, size=size)


def sci(x, d=0):
    """2.4e-04 -> '2 × 10^{−4}' (markup for a real superscript)."""
    m, e = ("%.*e" % (d, x)).split("e")
    return "%s × 10^{%s}" % (m, ("−" if int(e) < 0 else "") + str(abs(int(e))))


def f3(x, d=3):
    return ("%." + str(d) + "f") % x


def pm(v, d=3):
    return ("%." + str(d) + "f ± %." + str(d) + "f") % (v[0], v[1])


# ====================================================================== front matter
TITLE = "Damage–repair kinetics and phase stiffness in a site-diluted XY model: implications for irradiated superconductors"
p = doc.add_paragraph()
H.add_rich(p, TITLE, size=16, bold=True)
for line in ("Leon Sandler", "Independent researcher, Northbrook, Illinois, USA", "E-mail: sandler.leon@gmail.com",
             "ORCID: 0009-0007-4584-808X"):
    q = doc.add_paragraph()
    H.add_rich(q, line, size=10.5)
    q.paragraph_format.space_after = Pt(0)
doc.add_paragraph()


# ====================================================================== body, in document order
def run(name, part=None):
    code = open(os.path.join(HERE, name), encoding="utf-8").read()
    if part is not None:
        code = part(code)
    tmp = os.path.join(OUT, "_exec_" + name)             # so that tracebacks show the executed lines
    open(tmp, "w", encoding="utf-8").write(code)
    exec(compile(code, tmp, "exec"), globals())


MARK33 = "# ---------------------------------------------------------------------------------- 3.3"
run("part_abstract.py")
run("part_intro.py")
run("part_methods.py")
run("part_results_a.py", lambda c: c[:c.index(MARK33)])            # 3.1 (and the helpers used by 3.3)
run("part_results_b.py")                                            # 3.2
run("part_results_a.py", lambda c: c[:c.index("# ---------------------------------------------------------------------------------- 3.1")] + c[c.index(MARK33):])  # 3.3
run("part_results_e.py")                                            # 3.3 continued: size dependence and threshold
run("part_results_c.py")                                            # 3.4
run("part_results_d.py")                                            # 3.5
run("part_sec4.py")                                                 # 4
run("part_back.py")                                                 # 5, 6, declarations, references

doc.core_properties.author = "Leon Sandler"
doc.core_properties.title = TITLE
VERSION = os.environ.get("MS_VERSION", "2")
outp = os.path.join(OUT, "Damage_Repair_Phase_Stiffness_XY_v%s.docx" % VERSION)
doc.save(outp)
json.dump(TLAB_NEW, open(LABFILE, "w"))
print("saved", outp, "| figures:", FIGN[0], "| tables:", TABN[0], "| references:", len(CITE), "| labels:", TLAB_NEW,
      "(labels changed since the previous pass: %s)" % (TLAB_NEW != TLAB))
