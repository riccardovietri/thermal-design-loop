# Literature survey

## Second pass (2026-11-11, full-text sources)
Read in full text unless marked.

- **OCP OAI System Liquid Cooling Guidelines** (Meta, Microsoft, Intel, CoolerMaster, Wiwynn, Boyd and others, March 2023):
  1.5 LPM/kW as the typical target (10 K rise), range 1.25–2.0 LPM/kW; under 0.02 °C/W case-to-inlet with interface material for
  a 1 kW module on a 30 × 30 mm heater, PG25, 2–5 LPM, "ideal conditions"; test-vehicle lid recommended below 75 °C; cold-plate
  structures "could be straight fin, pin fin, and offset fin", internal features 0.1–1.0 mm.
  https://www.opencompute.org/documents/oai-system-liquid-cooling-guidelines-in-ocp-template-mar-3-2023-update-pdf
- **Ortega, Caceres, Uras, Kisitu, Chowdhury, Radmard, Heydari**, "Determination of the Thermal Performance Limits for Single Phase
  Liquid Cooling Using an Improved Effectiveness-NTU Cold Plate Model", ASME InterPACK2022-97421: skived copper plate, fin 0.2 mm,
  channel 0.2 mm, height 4 mm, length 43 mm, PG25 at 32 °C inlet, 1 kW; minimum measured 0.012 °C/W at 2.6 LPM (base to inlet);
  6.3 psi at 2.6 LPM for the whole plate; fully laminar (Re < 200). Plate width and base thickness are not given in the text.
  https://par.nsf.gov/servlets/purl/10435523
- **Ramakrishnan et al.** (SUNY Binghamton), cold plates for multi-chip server modules: 6.2–6.7 kPa at 1 LPM, DI water, whole plate
  including manifolds. https://par.nsf.gov/servlets/purl/10058150
- **CoolIT 4 kW plate** (vendor page): under 0.009 °C/W at 1.5 LPM per kW; 8 psi across the full loop including fittings and
  quick-disconnects; inlet temperature and footprint not stated.
- **OCP ACS Cold Plate Requirements** and **Cold Plate Cooling Loop Requirements Rev 2**: list what to publish (pressure-flow and
  resistance-flow curves); no numeric pressure-drop or resistance limit.
- No public case-temperature limit was found for H100 SXM, B200 or MI300-class parts. Package powers found: H100 SXM up to 700 W
  (vendor page), HGX B200 1000 W (Lenovo product guide), MI300X 750 W and MI355X 1400 W (AMD datasheets).
- **Offset strip fins:** the Manglik & Bergles (1995) correlations (constants confirmed from two secondary full texts; original not
  read) are fitted mostly to air data at moderate to high Reynolds number; Vangeffelen, Buckinx, Vetrano, Baelmans (KU Leuven,
  arXiv:2107.10719 and arXiv:2204.09605) report that in steady laminar micro-channels they under-predict friction by 20% on average
  and Nusselt number by 70% on average, and give replacement correlations of the form f = c0/Re + c1 and Nu = d0 + d1·Re (reference
  length is the strip length; read their definitions before use).
- **Wavy channels:** Kundu (2001, MS thesis, University of Cincinnati, supervisor Manglik), 2D and 3D sinusoidal corrugated ducts at
  Pr = 1: in 2D, friction and Nusselt number stay essentially flat at Re 50–100 and rise only from about Re 200; 3D adds
  longitudinal vortices from the side walls. Ghorbani et al. (Scientific Reports 12:5882, 2022), 3D conjugate simulations at
  Re 300–900: Nusselt number up to 3.95× and friction 3.15× straight at Re 500.
- No source was found comparing these planforms with straight fins at equal pumping power or equal pressure drop at Re 100–200.
- **Manufacturing:** Boyd fabrication guide: copper skived fins 0.1–0.8 mm thick, 0.1 mm minimum gap. No source gives a volume
  process or feature limit for offset-strip, wavy or pin-fin copper at 0.1–0.3 mm.
- Not found: a source for 0.03 K/W as a power-electronics cold-plate minimum (nearest: under 0.04 °C/W for data-centre liquid
  cooling, Amkor IMAPS slides 2021).

## First pass (2026-11-10, automated web survey)

Produced by an automated web survey (about 20 searches, 20 page fetches). **Confidence is uneven.** Many PDFs could not be opened
(OCP specifications, the Stanford and Kong manifold papers, the Boyd guide), so entries marked *snippet* come from search-result
excerpts and have not been checked against the full text. Do not cite a *snippet* number externally without reading the source.

### What it means for this project
- **Pressure budget.** No OCP or ASHRAE document that could be read gives a per-cold-plate limit. The only per-plate figure is a
  vendor listing: under 20 kPa at about 1.8 L/min for a GB200-class plate. The 6 kPa budget therefore remains an allocation. (superseded by second pass)
- **Equal-pressure comparisons are missing from the literature found.** No study gave numbers for any topology against straight
  microchannels at equal pressure drop or pumping power. This study's Pareto comparison addresses that directly.
- **Jets vs microchannels depends on size.** Lee & Vafai (1999): microchannels preferred below about 70 × 70 mm, jets comparable
  or better above (abstract only). Robinson (2009): jets run at lower pressure and higher flow, microchannels the reverse.
  Both agree with the slot-jet sweep here.
- **Manifold microchannels** are the best-documented high-performance option, but on heaters of about 1 cm² or less (*snippet*).
- **Reduced-order models for manifold microchannels exist**, which supports modelling them in 2D here (table below).

### Topologies
| Topology | Reported result | Pressure / flow | Basis | Source |
|---|---|---|---|---|
| Manifold microchannel, CFD design, water | 700 W/cm² at 0.023 K/W | 22 kPa, 1.3 L/min | Simulated, *snippet* | Kong et al., J. Electron. Packag. 2025 |
| Manifold microchannel, copper inverse-opal | 0.068 cm²·K/W convective, above 1 kW/cm² | 32 kPa total | Measured, *snippet* | Kong et al. 2024 |
| Manifold microchannel | 15.5 W/K per cm² | 0.005 W pumping power | Measured, as cited in a review | Copeland 2005 (citing Muller/Frechette) |
| Jet array vs microchannel | Both remove 250 W/cm² at 85 °C under 0.1 W pumping | Jets: lower ΔP, higher flow | Model, abstract | Robinson, IEEE Trans. CPMT 2009 |
| Jet interlayer, 1 cm² | 180 W/cm² at 100 µm pitch | 1 bar, 300 mL/min | Simulated, *snippet* | Brunschwiler et al. 2009 |
| Microchannel + pin-fin hybrid | 0.208 K·cm²/W | about 10 kPa | *snippet*, source not confirmed | — |
| Wavy channels | Nu +32% vs straight at 4 m/s | Penalty smaller than gain (qualitative) | Simulated, *snippet* | — |
| Gyroid lattice | About 2.4× lower resistance than pin fins | ΔP +63% for 10 → 6 mm cells | Mixed, *snippet* | — |

### Manufacturing limits
| Process | Limit reported | Source |
|---|---|---|
| Skiving, copper | Fin 0.1–0.8 mm, minimum gap 0.1 mm | Boyd fabrication guide (*snippet*) (superseded by second pass) |
| Skiving, other vendors | Fin 0.15–1.0 mm, pitch 0.4–3.5 mm; aspect ratio 12–15 | Vendor catalogues, a patent; values disagree |
| CNC milling | Channel and fin 0.5 mm | Copeland 2005 |
| Copper mini-channel optimum | Fin and gap 100–200 µm, height 1–2 mm | Ortega & Caceres, ASME 2022 (*snippet*) (superseded by second pass) |
| Etching + diffusion bonding, laser cutting, additive | Not found | — |

### System-level requirements
- ASHRAE liquid-cooling classes W17, W27, W32, W40, W45, W+; a 40 °C supply sits in W40 (vendor summary, not ASHRAE's text).
- OCP Cold Plate Test Specification V1.0.0 defines pressure-drop and flow test methods, apparently not limits (*snippet*). (superseded by second pass)
- NVIDIA GB200 NVL72, secondary sources that conflict: inlet 18–32 °C, 20–25 °C or up to 45 °C; rack flow 80 to 130 L/min; 1.5 bar system pressure.
- Not found: a primary-source flow-per-kW figure; case or junction temperature limits (not searched). (superseded by second pass)

### Reduced-order models
| Model | Agreement | Checked against | Source |
|---|---|---|---|
| Manifold microchannel correlations, fit to 490 CFD cases | Friction within 0.7%; Nu mean 0.4%, worst 11–12% | CFD only | Yuruker, Mandel, Shooshtari, Ohadi, ITherm 2019 / Electronics Cooling 2022 |
| Porous-zone manifold model | Pressure drop within 5% | Experiment; thermal side not validated (*snippet*) | 2024 conference paper |
| Porous-medium microchannel layer | Not retrieved | 3D single-channel CFD (*snippet*) | Nonino & Savino 2025 |

### Sources read in full or in part
- Copeland, "Review of Low Profile Cold Plate Technology for High Density Servers", Electronics Cooling, 2005: https://www.electronics-cooling.com/2005/05/review-of-low-profile-cold-plate-technology-for-high-density-servers/
- Yuruker et al., "Heat Transfer and Pressure Drop Correlations for Manifold Microchannel Heatsinks", Electronics Cooling, 2022: https://www.electronics-cooling.com/2022/10/heat-transfer-and-pressure-drop-correlations-for-manifold-microchannel-heat-sinks/
- Robinson, IEEE Trans. CPMT 32(2), 2009 (abstract): https://edepositireland.ie/handle/2262/31928
- Lee & Vafai, Int. J. Heat Mass Transfer 42(9), 1999 (abstract): https://portal.mardi4nfdi.de/wiki/Item:Q1973216
- Fibermall GB200 cold plate listing (vendor): https://www.fibermall.com/blog/nvidia-gb200-liquid-cooling-plate.htm
- QCT Qoolrack white paper for GB200 NVL72 (vendor): https://blog.qct.io/wp-content/uploads/2025/04/QCT-Qoolrack-Stand-Alone_Advanced-Liquid-Cooling-for-NVIDIA-GB200-NVL72-Systems.pdf
- Kong et al., J. Electron. Packag. 2025 (snippet): https://atsla.cau.ac.kr/wp-content/uploads/2025/09/73-Kong-et-al-JEP2025.pdf
- OCP Cold Plate Test Specification V1.0.0 (unreadable PDF): https://www.opencompute.org/wp-content/uploads/2026/10/OCP-Cold-Plate-Test-Specification-V1.0.0.pdf

### Follow-up
Read the full text of Kong 2024/2025, Jung 2020 and Robinson 2009; obtain the OCP documents in readable form; find a primary
NVIDIA or partner coolant specification; search for accelerator case-temperature limits and for equal-pumping-power comparisons.
(OCP documents, case-temperature limits and equal-pumping-power comparisons: superseded by second pass.)

## Third pass (2026-11-15, entrance effects, issue #12)
Full text read unless marked. PDFs fetched from the publishers' or authors' institutional pages and read with pypdf.

- **Muzychka, Y.S. and Yovanovich, M.M. (2002)**, "Laminar flow friction and heat transfer in non-circular ducts and channels,
  Part I: hydrodynamic problem", in *Compact Heat Exchangers: A Festschrift on the 60th Birthday of Ramesh K. Shah*, Grenoble,
  August 2002, pp. 123–130. Full text: Memorial University repository, https://hdl.handle.net/20.500.14783/9794 (copy of
  http://www.mhtlab.uwaterloo.ca/pdf_papers/mhtl02-11.pdf). Used: apparent friction f_app·Re_√A = [(f·Re_√A)² + (3.44/√L+)²]^½
  (Eq. 43), L+ = L/(√A·Re_√A), f·Re_√A from the first term of the rectangular-duct series (Eq. 32); Table 1 (exact and Eq. 32
  values for 0.001 ≤ ε ≤ 1); worked values L+_h = 0.059 (ε = 1) and 0.00083 (ε = 0.01) from Eq. 44. Stated accuracy ±10% or
  better for most duct shapes; parallel plates modelled as ε = 0.01. The model has no separate K(∞) term: Siegel's short-duct
  solution already carries the momentum change of the accelerating core. Source of Shah & London data: their book, not read here.
- **Muzychka, Y.S. and Yovanovich, M.M. (2002)**, "... Part II: thermal problem", same volume, pp. 131–139. Full text: Memorial
  University repository, https://research.library.mun.ca/2486 (now
  https://memorial.scholaris.ca/items/6a2d925d-d1d4-4904-b6a2-ef3e49d4c888). Used: Graetz model (Eq. 50) and combined
  entry model (Eqs. 51–52, m = 2.27 + 1.65 Pr^(1/3)), developed limit Nu_√A = C1·f·Re_√A/(8√π·ε^γ) (Eq. 31), flat-plate function
  f(Pr) for uniform flux (Eq. 48), constants of Table 2 (uniform flux: C1 = 3.66, C3 = 0.501; local C2 = C4 = 1, mean 3/2 and 2;
  γ = 1/10 upper bound, −3/10 lower bound); worked values z*_t = 0.0366 (H) and 0.0358 (T) for the circular duct (Eq. 53).
  Stated accuracy ±12% (Graetz) and ±15% (combined entry); valid 0.1 < Pr < ∞.
- **dos Reis, E. and Altemani, C.A.C. (1999)**, "Pressure drop and heat transfer in airflow cooled finned heat sinks using
  correlations", COBEM 1999 (15th Brazilian Congress of Mechanical Engineering). Full text:
  https://abcm.org.br/anais/cobem/1999/english/..\pdf\AACGEC.pdf (the URL as indexed). Used: polynomial fits to the laminar abrupt-contraction and
  abrupt-expansion coefficients of Kays (1950) for a flat duct and a square duct, interpolated linearly in aspect ratio (Eq. 5,
  Table 1); the pressure-drop sum ΔP = ρU²/2 (Kc + 4f·L/Dh + Ke) with an apparent friction factor (Eq. 2). Secondary source: the
  original Kays (1950), "Loss coefficients for abrupt changes in flow cross section with low Reynolds number flow in single and
  multiple-tube systems", Trans. ASME 72, 1067–1074, was **not** obtained, nor Kays & London, *Compact Heat Exchangers* (Internet
  Archive copy is lending-only). The subscripts of the Table 1 coefficients are lost in the PDF text; the order A0, A1, A2 is
  inferred from the row order and checked: the flat-duct expansion fit equals the momentum balance with a parabolic exit profile.
- **Qu, W. and Mudawar, I. (2002)**, "Experimental and numerical study of pressure drop and heat transfer in a single-phase
  micro-channel heat sink", Int. J. Heat Mass Transfer 45, 2549–2565, https://engineering.purdue.edu/mudawar/files/articles-all/2002/2002_1.pdf.
  Full text. Includes contraction and expansion losses in a measured micro-channel plate (Eqs. 18–22) with coefficients
  "close to unity" for their geometry; no correlation given. Not used for numbers.
- **Phillips, R.J. (1988)**, "Microchannel heat sinks", Lincoln Laboratory Journal 1(1),
  https://archive.ll.mit.edu/publications/journal/pdf/vol01_no1/1.1.3.microchannelheatsinks.pdf. Full text. Laminar apparent
  friction and developing Nu for rectangular channels are given only as figures; loss coefficients set to zero in the examples.
  Not used for numbers.
- **Johnson, B.J. (2004)**, "Pressure losses for fluid flow through abrupt area contraction in compact heat exchangers",
  undergraduate research report, University of Washington, https://faculty.washington.edu/finlayso/che499/Bryan_Sp04.pdf. Full
  text. 2D FEMLAB contraction coefficients against the Kays & London "Laminar" line: 17–37% different at Re = 100 depending on σ.
  Unreviewed student work; noted as a warning that the laminar Kc depends on Re, not used for numbers.
- **Wanted, not obtained:** Shah & London (1978) tables of apparent friction and developing Nu for parallel plates and
  rectangular ducts (book); Shah (1978), "A correlation for laminar hydrodynamic entry length solutions for circular and non-circular ducts", J. Fluids
  Eng. 100, 177–179; Kays (1950);
  Muzychka & Yovanovich (2009), J. Fluids Eng. 131, 111105 (no open copy found). Phillips's MIT thesis (1987, DSpace handle 1721.1/14921)
  could not be fetched (rate-limited).
