# -*- coding: utf-8 -*-
# Section 1 (exec'd by build_manuscript.py)
HD("1. Introduction and superconducting motivation")
P("High-temperature superconducting (HTS) magnets are central to compact fusion concepts [[sorbom2015,creely2020]]. In service the conductor is exposed "
  "to fast neutrons, which displace atoms and degrade its superconducting properties [[zinkle2014]]. For REBCO coated conductors, fast-neutron "
  "irradiation first raises the critical current, which passes through a temperature-dependent maximum and then declines [[prokopec2015,fischer2018]]; "
  "at 50 K it falls below its pristine value at fluences of the order of 10^{22} m^{−2}, and the *n*-value decreases even at low "
  "fluence [[fischer2019,unterrainer2022]]. The transition temperature, by contrast, decreases linearly with the calculated number of displacements per atom (dpa), "
  "at the same rate for neutrons and for helium ions [[adams2023]]. Annealing partly reverses the damage: about 25% of the decrease of *T*_{c} "
  "was recovered at 150 °C and about 60% at 400 °C, and the degradation of *T*_{c} and of the critical current density were found to be closely "
  "related, “likely due to the expected loss of superfluid density” [[unterrainer2022]]. The engineering question is therefore not only how "
  "fast damage accumulates but whether superconducting order persists when damage and recovery act together.")
P("The superfluid density is the physical counterpart of the phase stiffness of the order parameter. In two-dimensional superconductors the "
  "stiffness controls the Berezinskii–Kosterlitz–Thouless (BKT) transition [[kt1973,nk1977,beasley1979,halperin1979,hebard1980]], and in any "
  "superconductor it sets the London penetration depth [[tinkham2004]]. A minimal model of how damage and recovery compete for phase coherence is "
  "therefore a two-dimensional XY model of the local phase, in which sites are damaged (pair breaking) and repaired. Diluted XY models "
  "and their BKT transitions have been studied [[costa2014]]; what is new here is the kinetics. In the model used, damage is faster where phase order "
  "is already disrupted, as would be the case if disorder lowered the pair amplitude locally, and repair acts at a stochastic rate. The model is "
  "implemented in the open-source MyUncle framework [[sandler_myuncle]].")
P("We ask four questions. Does state-dependent damage make the response to the repair rate bistable and hysteretic, as a mean-field reduction "
  "suggests? Does recovery restore genuine phase *stiffness*, not merely alignment in an external field? Is the repair-to-damage ratio ρ the "
  "right control parameter, or does the competition between the phase-relaxation and damage–repair timescales matter? And what can such a model say "
  "about the experiments above, given that it contains no microscopic pairing physics, no vortex pinning and no transport current, and has "
  "not been calibrated to any conductor? We answer each with the numerical evidence stated together with its limits.")
