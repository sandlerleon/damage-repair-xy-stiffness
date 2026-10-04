# -*- coding: utf-8 -*-
"""Reference list. Every entry with a DOI was resolved against Crossref (refs_cache.json holds the records; run with
--refresh to re-harvest) and is formatted from the cached record; books and non-Crossref items are entered by hand."""
import io
import json
import os
import re
import sys
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "refs_cache.json")

DOI = {
    "sorbom2015": "10.1016/j.fusengdes.2015.07.008",
    "creely2020": "10.1017/S0022377820001257",
    "zinkle2014": "10.1146/annurev-matsci-070813-113627",
    "prokopec2015": "10.1088/0953-2048/28/1/014005",
    "fischer2018": "10.1088/1361-6668/aaadf2",
    "adams2023": "10.1088/1361-6668/aced9e",
    "unterrainer2022": "10.1088/1361-6668/ac4636",
    "beasley1979": "10.1103/PhysRevLett.42.1165",
    "hebard1980": "10.1103/PhysRevLett.44.291",
    "halperin1979": "10.1007/BF00116988",
    "kt1973": "10.1088/0022-3719/6/7/010",
    "nk1977": "10.1103/PhysRevLett.39.1201",
    "wm1988": "10.1103/PhysRevB.37.5986",
    "hasenbusch2005": "10.1088/0305-4470/38/26/003",
    "newman2000": "10.1103/PhysRevLett.85.4104",
    "costa2014": "10.1088/1742-6596/487/1/012008",
    "rao1990": "10.1103/PhysRevB.42.856",
    "efron1979": "10.1214/aos/1176344552",
    "madras1988": "10.1007/BF01022990",
}
TITLE = {   # sentence-case titles, and records whose Crossref title loses its formulas
    "sorbom2015": "ARC: a compact, high-field, fusion nuclear science facility and demonstration power plant with demountable magnets",
    "zinkle2014": "Designing radiation resistance in materials for fusion energy",
    "unterrainer2022": "Recovering the performance of irradiated high-temperature superconductors for use in fusion magnets",
    "adams2023": "Comparing neutron and helium ion irradiation damage of REBa2Cu3O7-δ coated conductor using x-ray absorption spectroscopy",
    "beasley1979": "Possibility of vortex-antivortex pair dissociation in two-dimensional superconductors",
    "hebard1980": "Evidence for the Kosterlitz-Thouless transition in thin superconducting aluminum films",
    "halperin1979": "Resistive transition in superconducting films",
    "kt1973": "Ordering, metastability and phase transitions in two-dimensional systems",
    "nk1977": "Universal jump in the superfluid density of two-dimensional superfluids",
    "wm1988": "Monte Carlo determination of the critical temperature for the two-dimensional XY model",
    "hasenbusch2005": "The two-dimensional XY model at the transition temperature: a high-precision Monte Carlo study",
    "newman2000": "Efficient Monte Carlo algorithm and high-precision results for percolation",
    "costa2014": "Kosterlitz-Thouless transition: the diluted XY model",
    "efron1979": "Bootstrap methods: another look at the jackknife",
    "madras1988": "The pivot algorithm: a highly efficient Monte Carlo method for the self-avoiding walk",
    "rao1990": "Magnetic hysteresis in two model spin systems",
}
JOURNAL = {
    "Fusion Engineering and Design": "Fusion Eng. Des.", "Journal of Plasma Physics": "J. Plasma Phys.",
    "Annual Review of Materials Research": "Annu. Rev. Mater. Res.", "Superconductor Science and Technology": "Supercond. Sci. Technol.",
    "Physical Review Letters": "Phys. Rev. Lett.", "Physical Review B": "Phys. Rev. B", "Journal of Low Temperature Physics": "J. Low Temp. Phys.",
    "Journal of Physics C: Solid State Physics": "J. Phys. C: Solid State Phys.", "Journal of Physics A: Mathematical and General": "J. Phys. A: Math. Gen.",
    "Journal of Physics: Conference Series": "J. Phys.: Conf. Ser.", "The Annals of Statistics": "Ann. Stat.",
    "Journal of Statistical Physics": "J. Stat. Phys.",
}
PAGES = {"efron1979": "1-26", "creely2020": "865860502"}       # fields Crossref leaves empty or that need the article number
MANUAL = {
    "fischer2019": "Fischer D X 2019 Effect of neutron radiation damage on coated conductors for fusion magnets PhD Thesis TU Wien "
                   "https://doi.org/10.34726/hss.2019.27911",
    "tinkham2004": "Tinkham M 2004 Introduction to Superconductivity 2nd edn (Mineola, NY: Dover)",
    "sandler_myuncle": "Sandler L 2026 MyUncle: a compact maintenance-lattice framework for order-persistence simulations across disciplines "
                       "(Zenodo) https://doi.org/10.5281/zenodo.21223569",
    "sandler_physa": "Sandler L 2026 Self-maintained order and hysteretic collapse in a non-equilibrium rotational lattice (Zenodo preprint) "
                     "https://doi.org/10.5281/zenodo.21210708",
}


def refresh():
    sys.path.insert(0, r"C:\YouTube\tncc-framework\references")
    import harvest as h
    db = json.load(io.open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}
    for k, d in DOI.items():
        if k in db and "--refresh" not in sys.argv:
            continue
        it = h.by_doi(d)
        if not it:
            raise SystemExit("Crossref has no record for %s" % d)
        db[k] = h.record(it)
    io.open(CACHE, "w", encoding="utf-8").write(json.dumps(db, indent=1, ensure_ascii=False))
    return db


def load():
    if not os.path.exists(CACHE):
        return refresh()
    return json.load(io.open(CACHE, encoding="utf-8"))


def initials(given):
    out = []
    for part in re.split(r"[\s]+", (given or "").strip()):
        if part:
            out.append("".join(s[0].upper() for s in part.split("-") if s))
    return " ".join(out)


def entry(key, db):
    """IOP numerical style: Surname A B, Surname A B et al Year Title Journal Volume Pages DOI."""
    if key in MANUAL:
        return MANUAL[key]
    r = db[key]
    au = ["%s %s" % (a["family"].replace("\ufffd", "ö"), initials(a["given"])) for a in r["authors"]]
    au = ", ".join(au) if len(au) <= 6 else ", ".join(au[:3]) + " et al"
    jn = {unicodedata.normalize("NFC", k): v for k, v in JOURNAL.items()}.get(unicodedata.normalize("NFC", r["container"]), r["container"])
    pages = PAGES.get(key) or (r.get("page") or "")
    pages = pages.replace("-", "\u2013")
    vol = r.get("volume") or ""
    assert vol and pages, "incomplete reference %s: volume=%r pages=%r" % (key, vol, pages)
    return "%s %s %s %s %s %s %s https://doi.org/%s" % (au, r["year"], TITLE.get(key, r["title"]) + ".", jn, vol, pages, "", DOI[key]) \
        if False else "%s %s %s %s %s, %s. https://doi.org/%s" % (au, r["year"], TITLE.get(key, r["title"]) + ".", jn, vol, pages, DOI[key])


if __name__ == "__main__":
    db = refresh()
    for k in DOI:
        print(k, "|", entry(k, db))
