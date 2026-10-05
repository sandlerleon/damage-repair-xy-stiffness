# -*- coding: utf-8 -*-
"""Publish the two Zenodo records reserved by zenodo_reserve.py:

  software     the tagged GitHub release as a zip (MIT)
  preprint     the manuscript (docx and pdf) (CC BY 4.0)

The DOIs were reserved first so that they could be written into the manuscript. Headline numbers in the descriptions are read
from results/results.json. The token is read from ZENODO_TOKEN and never written to disk.

    python zenodo_publish.py software|preprint [--version=1.0.0] [--ms=1] [--dry]
"""
import json
import os
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request

TOKEN = os.environ.get("ZENODO_TOKEN")
if not TOKEN:
    raise SystemExit("ZENODO_TOKEN is not set in the environment")
API = "https://zenodo.org/api"
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
STATE = os.path.join("C:" + os.sep, "YouTube", "_dxy_zenodo_state.json")
VERSION = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--version=")), "1.0.0")
MS = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--ms=")), "1")
TAG = "v" + VERSION
DRY = "--dry" in sys.argv
GITHUB = "https://github.com/sandlerleon/damage-repair-xy-stiffness"
CREATORS = [{"name": "Sandler, Leon", "affiliation": "Independent Researcher", "orcid": "0009-0007-4584-808X"}]
TITLE_PAPER = "Damage-repair kinetics and phase stiffness in a site-diluted XY model: implications for irradiated superconductors"
TITLE_CODE = ("Damage-repair kinetics and phase stiffness in a site-diluted XY model: kinetic Monte Carlo engine, twist-response "
              "stiffness, analysis and manuscript")
KEYWORDS = ["superconductivity", "kinetic Monte Carlo", "XY model", "phase stiffness", "helicity modulus", "twist response",
            "irradiation damage", "annealing", "REBCO", "fusion magnets", "Berezinskii-Kosterlitz-Thouless", "site dilution"]
R = json.load(open(os.path.join(REPO, "results", "results.json"), encoding="utf-8"))

ABOUT = """<p><strong>A computational study of a phenomenological model. No experiment was performed and the model is not calibrated to any
conductor.</strong> Prepared for submission to <em>Superconductor Science and Technology</em>. A two-dimensional XY lattice of the local
superconducting phase, in which sites are damaged (pair breaking) at a rate that rises with local phase disorder and repaired at a stochastic rate,
is simulated by kinetic Monte Carlo (engine verified bit for bit against the MyUncle framework,
<a href="https://doi.org/10.5281/zenodo.21223569">10.5281/zenodo.21223569</a>).</p>
<p><strong>What the study contains.</strong> (i) The mean-field reduction re-derived from the implemented update for the damaged fraction, with
the bistable windows of both the rate equation and the discrete-time map, and a test of its closure against the lattice. (ii) Two-start tests and
quasi-static ramps at baseline field, with all interval exceptions listed and Holm-corrected. (iii) A twist-response measurement of the phase
stiffness of the non-equilibrium steady state (with winding numbers recorded), compared with the equilibrium fluctuation formula, and phase
correlations C(r), for repair rules that preserve the U(1) symmetry. (iv) The competition between damage rate and phase relaxation
(timescale parameter) and the discrete-time effects. (v) Finite-size scaling of the quenched dilution transition to L = 128, clustered damage,
the BKT line T_BKT(f), and an illustrative mapping to displacement damage and annealing.</p>"""
DESC_CODE = ABOUT + """<p>Contents: the engine, its test suite, the numerical studies (<code>code/study.py</code>, including the twist-response validation), the analysis and figure scripts,
the manuscript build script, all raw run outputs, the figures and the manuscript. Every run is seeded from SeedSequence([master, code]) and
can be repeated alone. Manuscript preprint: <a href="https://doi.org/{PP}">{PP}</a>.</p>"""
DESC_PAPER = ABOUT + """<p>Code, raw data and analysis: <a href="%s">%s</a>, archived at
<a href="https://doi.org/{SW}">{SW}</a>.</p>""" % (GITHUB, GITHUB)


def req(method, url, data=None, headers=None, raw=None):
    h = {"Authorization": "Bearer " + TOKEN}
    if headers:
        h.update(headers)
    body = raw if raw is not None else (json.dumps(data).encode() if data is not None else None)
    if data is not None and raw is None:
        h["Content-Type"] = "application/json"
    r = urllib.request.Request(url, data=body, headers=h, method=method)
    try:
        with urllib.request.urlopen(r, timeout=600) as resp:
            t = resp.read()
            return json.loads(t) if t else {}
    except urllib.error.HTTPError as e:
        raise SystemExit("%s %s -> %s\n%s" % (method, url, e.code, e.read().decode()[:800]))


def upload(bucket, path, name):
    with open(path, "rb") as fh:
        req("PUT", "%s/%s" % (bucket, urllib.parse.quote(name)), raw=fh.read(), headers={"Content-Type": "application/octet-stream"})
    print("   uploaded %-62s %9.1f kB" % (name, os.path.getsize(path) / 1024.0))


def finish(did, meta):
    req("PUT", "%s/deposit/depositions/%s" % (API, did), data={"metadata": meta})
    print("   metadata written")
    if DRY:
        print("   DRY RUN - draft %s left unpublished" % did)
        return
    pub = req("POST", "%s/deposit/depositions/%s/actions/publish" % (API, did))
    rec = req("GET", "%s/records/%s" % (API, pub["id"]))
    print("   PUBLISHED  DOI %s  concept %s" % (rec.get("doi"), rec.get("conceptdoi")))


def clear_inherited(d):
    """A new-version draft starts with the files of the previous version; remove them so that only this version's files remain."""
    for fid in d.get("inherited_files", []):
        req("DELETE", "%s/deposit/depositions/%s/files/%s" % (API, d["id"], fid))


NEWVER = ("<p><strong>Version %s.</strong> Revised after review: adds control simulations of the twist response (clean lattices, imposed windings, an independent "
          "static code; Section 2.5), qualifies the threshold brackets as finite-size estimates, presents the slow-damage convergence as numerical evidence, removes the "
          "interpretation of the sites-per-dpa factor, and qualifies the annealing analysis.</p>")


def software():
    st = json.load(open(STATE))
    d = st["software_" + VERSION] if "software_" + VERSION in st else st["software"]
    tmp = os.path.join(os.environ.get("TEMP", "."), "damage-repair-xy-stiffness-%s.zip" % VERSION)
    subprocess.check_call(["git", "-C", REPO, "archive", "--format=zip", "--prefix=damage-repair-xy-stiffness-%s/" % VERSION, "-o", tmp, TAG])
    print("=== software draft %s (reserved DOI %s)" % (d["id"], d["doi"]))
    clear_inherited(d)
    upload(d["bucket"], tmp, os.path.basename(tmp))
    meta = {"title": TITLE_CODE, "upload_type": "software", "description": (NEWVER % VERSION if VERSION != "1.0.0" else "") + DESC_CODE.replace("{PP}", st.get("publication_v" + MS, st["publication"])["doi"]),
            "creators": CREATORS, "keywords": KEYWORDS, "access_right": "open", "license": "mit-license", "version": VERSION, "language": "eng",
            "prereserve_doi": {"doi": d["doi"]},
            "related_identifiers": [{"identifier": GITHUB + "/tree/" + TAG, "relation": "isSupplementTo", "scheme": "url"},
                                    {"identifier": st.get("publication_v" + MS, st["publication"])["doi"], "relation": "isSupplementTo", "scheme": "doi"},
                                    {"identifier": "10.5281/zenodo.21223569", "relation": "isDerivedFrom", "scheme": "doi"}]}
    finish(d["id"], meta)


def preprint():
    st = json.load(open(STATE))
    d = st["publication_v" + MS] if "publication_v" + MS in st else st["publication"]
    print("=== preprint draft %s (reserved DOI %s)" % (d["id"], d["doi"]))
    clear_inherited(d)
    for name in ("Damage_Repair_Phase_Stiffness_XY_v%s.docx" % MS, "Damage_Repair_Phase_Stiffness_XY_v%s.pdf" % MS):
        upload(d["bucket"], os.path.join(REPO, "manuscript", name), name)
    meta = {"title": TITLE_PAPER, "upload_type": "publication", "publication_type": "preprint",
            "description": (NEWVER % ("v" + MS) if MS != "1" else "") + DESC_PAPER.replace("{SW}", st.get("software_" + VERSION, st["software"])["doi"]), "creators": CREATORS, "keywords": KEYWORDS, "access_right": "open",
            "license": "cc-by-4.0", "version": MS, "language": "eng", "prereserve_doi": {"doi": d["doi"]},
            "related_identifiers": [{"identifier": st.get("software_" + VERSION, st["software"])["doi"], "relation": "isSupplementedBy", "scheme": "doi"},
                                    {"identifier": GITHUB, "relation": "isSupplementedBy", "scheme": "url"},
                                    {"identifier": "10.5281/zenodo.21210708", "relation": "references", "scheme": "doi"},
                                    {"identifier": "10.5281/zenodo.21223569", "relation": "references", "scheme": "doi"}]}
    finish(d["id"], meta)


if __name__ == "__main__":
    what = [a for a in sys.argv[1:] if not a.startswith("--")]
    if what == ["software"]:
        software()
    elif what == ["preprint"]:
        preprint()
    else:
        print(__doc__)
