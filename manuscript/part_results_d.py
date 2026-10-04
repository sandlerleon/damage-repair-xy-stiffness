# -*- coding: utf-8 -*-
# Section 3.5: static dilution, extrapolation, clustered damage (exec'd by build_manuscript.py)
FS_ = R["fss"]["0.35"]
CLU = R["cluster"]
PL_ = FS_["per_L"]
Lsz = [L for L in ("16", "24", "32", "48", "64", "96") if L in PL_]
EXA = FS_["extrap_all"]
SEN = FS_["sensitivity_lnL2"]
WMs = FS_["weber_minnhagen"]

HD("3.5 Static dilution, finite-size extrapolation and clustered damage", 2)
cis = [(PL_[L]["ci95"][1] - PL_[L]["ci95"][0]) / 2 for L in Lsz]
P("The static limit, in which a fraction *f* of the sites is removed at random and the phases equilibrate (*g* = 0 on removed sites, *T* = 0.35, "
  "equation (7), which is legitimate in equilibrium), provides the reference for the dynamic results and the input to the mapping of Section 4. "
  "Υ(*f*) was computed for *L* = 16 to 96 (and 128 near the crossing) with 6–12 disorder realizations per point, 2000 equilibration and 4000 "
  "measurement sweeps, and a grid of 0.01 near the crossing (Figure 7a). The crossing with the Nelson–Kosterlitz value 2*T*/π = %s [[nk1977]] "
  "is at *f*_{KT}(*L*) = %s for *L* = %s, each with a 95%% bootstrap interval over realizations of half-width %s–%s (@T:cross@). The scatter of the "
  "crossings over sizes, standard deviation %s, is smaller than the uncertainty of each one and must not be used as an error estimate; the "
  "uncertainty of the limit comes from the extrapolation. Extrapolating *f*_{KT}(*L*) linearly in 1/ln²*L* gives *f*_{KT}(∞) = %s (95%% interval %s–%s), "
  "and a Weber–Minnhagen fit [[wm1988]] of Υ_{L} = (2*T*/π)[1 + 1/(2 ln *L* + *C*)] to the sizes *L* ≥ 24 gives %s (%s–%s; χ² = %s for %d degrees of "
  "freedom). The estimate is stable against the smallest lattices used: excluding *L* = 16, then also 24 and 32, gives %s, %s and %s, with intervals that widen "
  "as fewer sizes remain (@T:cross@). We quote *f*_{KT} = %s ± 0.01 as a numerical estimate for *T* = 0.35, from sizes up to 96 and 6–12 realizations, "
  "not as a precisely established critical value. It lies well below the site-percolation threshold 1 − *p*_{c} = 0.407 [[newman2000]]: the "
  "intact sublattice still percolates, and the loss of stiffness is a phase-fluctuation effect, as in the BKT picture."
  % (f3(FS_["NK"]), ", ".join(f3(PL_[L]["f_KT"]) for L in Lsz), ", ".join(Lsz), f3(min(cis), 3), f3(max(cis), 3), f3(EXA["scatter_sd_over_sizes"], 4),
     f3(EXA["f_inf"]), f3(EXA["ci95"][0]), f3(EXA["ci95"][1]), f3(WMs["24"]["f_KT"]), f3(WMs["24"]["ci95"][0]), f3(WMs["24"]["ci95"][1]),
     f3(WMs["24"]["chi2"], 2), WMs["24"]["dof"], f3(SEN["24"]["f_inf"]), f3(SEN["32"]["f_inf"]), f3(SEN["48"]["f_inf"]), f3(round(EXA["f_inf"], 2), 2)))
rows = [["*L*", "*f*_{KT}(*L*)", "95% CI of the crossing", "realizations"]]
for L in Lsz:
    rows.append([L, f3(PL_[L]["f_KT"], 4), "%s–%s" % (f3(PL_[L]["ci95"][0], 3), f3(PL_[L]["ci95"][1], 3)), str(PL_[L]["n_real"])])
rows.append(["∞, 1/ln²*L*, all sizes", f3(EXA["f_inf"], 4), "%s–%s" % (f3(EXA["ci95"][0], 3), f3(EXA["ci95"][1], 3)), "sizes 16–96"])
for k in ("24", "32", "48"):
    rows.append(["∞, 1/ln²*L*, *L* ≥ %s" % k, f3(SEN[k]["f_inf"], 4), "%s–%s" % (f3(SEN[k]["ci95"][0], 3), f3(SEN[k]["ci95"][1], 3)), "sizes %s–96" % k])
rows.append(["∞, Weber–Minnhagen, *L* ≥ 24", f3(WMs["24"]["f_KT"], 4), "%s–%s" % (f3(WMs["24"]["ci95"][0], 3), f3(WMs["24"]["ci95"][1], 3)), "sizes 24–96"])
TAB(rows, "Crossing of the helicity modulus with 2*T*/π under random dilution (*T* = 0.35): the crossing at each size with its own bootstrap interval, "
          "and the extrapolations to *L* → ∞ with their sensitivity to the smallest sizes used.", widths=[2.3, 1.2, 1.9, 1.2], label="cross")
cl_f = {rad: CLU["%d|48" % rad]["f_KT"] for rad in (0, 1, 2, 3)}
cl_ci = {rad: CLU["%d|48" % rad]["ci95"] for rad in (0, 1, 2, 3)}
P("Irradiation does not remove independent sites: displacement cascades produce extended defects. We compared random dilution with compact "
  "clusters of 5, 13 and 29 removed sites (discs of radius 1, 2 and 3, placed at random positions up to the same removed fraction *f*; *L* = 32 and 48, 10 "
  "realizations, *T* = 0.35). At *L* = 48 the crossing moves from %s for random sites to %s (%s–%s), %s (%s–%s) and %s (%s–%s) for the three cluster sizes "
  "(Figure 7c), towards the percolation value 0.407, and the same ordering holds at *L* = 32 (%s, %s, %s, %s). At the same removed fraction, clustered damage "
  "therefore destroys the stiffness far less than independent point damage, and the mapping of Section 4.2, which assumes random point damage, "
  "gives the number of non-superconducting sites per dpa only for that assumption; it does not distinguish a cascade of many sites that each count little "
  "from fewer isolated ones."
  % (f3(cl_f[0], 3), f3(cl_f[1], 3), f3(cl_ci[1][0], 3), f3(cl_ci[1][1], 3), f3(cl_f[2], 3), f3(cl_ci[2][0], 3), f3(cl_ci[2][1], 3), f3(cl_f[3], 3),
     f3(cl_ci[3][0], 3), f3(cl_ci[3][1], 3), *[f3(CLU["%d|32" % r]["f_KT"], 3) for r in (0, 1, 2, 3)]))
FIG("Fig7_dilution.png",
    "Static dilution. (a) Helicity modulus Υ against the removed fraction *f* at *T* = 0.35 for *L* = 16–128 (95%% CI); dashed: 2*T*/π. (b) The "
    "crossing at each size with its bootstrap interval against 1/ln²*L*, the linear extrapolation (line, square), its sensitivity to dropping the "
    "smallest sizes (triangles, *L* ≥ 24, 32, 48) and, as a grey band, the scatter of the crossings over sizes, which is not an error estimate. (c) Υ(*f*) at *L* = 48 "
    "for random removal and for compact clusters of radius 1, 2 and 3 (5, 13 and 29 sites).")
