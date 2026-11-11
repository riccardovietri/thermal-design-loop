# Literature survey (2026-11-10)

Produced by an automated web survey (about 20 searches, 20 page fetches). **Confidence is uneven.** Many PDFs could not be opened
(OCP specifications, the Stanford and Kong manifold papers, the Boyd guide), so entries marked *snippet* come from search-result
excerpts and have not been checked against the full text. Do not cite a *snippet* number externally without reading the source.

## What it means for this project
- **Pressure budget.** No OCP or ASHRAE document that could be read gives a per-cold-plate limit. The only per-plate figure is a
  vendor listing: under 20 kPa at about 1.8 L/min for a GB200-class plate. The 6 kPa budget therefore remains an allocation.
- **Equal-pressure comparisons are missing from the literature found.** No study gave numbers for any topology against straight
  microchannels at equal pressure drop or pumping power. This study's Pareto comparison addresses that directly.
- **Jets vs microchannels depends on size.** Lee & Vafai (1999): microchannels preferred below about 70 × 70 mm, jets comparable
  or better above (abstract only). Robinson (2009): jets run at lower pressure and higher flow, microchannels the reverse.
  Both agree with the slot-jet sweep here.
- **Manifold microchannels** are the best-documented high-performance option, but on heaters of about 1 cm² or less (*snippet*).
- **Reduced-order models for manifold microchannels exist**, which supports modelling them in 2D here (table below).

## Topologies
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

## Manufacturing limits
| Process | Limit reported | Source |
|---|---|---|
| Skiving, copper | Fin 0.1–0.8 mm, minimum gap 0.1 mm | Boyd fabrication guide (*snippet*) |
| Skiving, other vendors | Fin 0.15–1.0 mm, pitch 0.4–3.5 mm; aspect ratio 12–15 | Vendor catalogues, a patent; values disagree |
| CNC milling | Channel and fin 0.5 mm | Copeland 2005 |
| Copper mini-channel optimum | Fin and gap 100–200 µm, height 1–2 mm | Ortega & Caceres, ASME 2022 (*snippet*) |
| Etching + diffusion bonding, laser cutting, additive | Not found | — |

## System-level requirements
- ASHRAE liquid-cooling classes W17, W27, W32, W40, W45, W+; a 40 °C supply sits in W40 (vendor summary, not ASHRAE's text).
- OCP Cold Plate Test Specification V1.0.0 defines pressure-drop and flow test methods, apparently not limits (*snippet*).
- NVIDIA GB200 NVL72, secondary sources that conflict: inlet 18–32 °C, 20–25 °C or up to 45 °C; rack flow 80 to 130 L/min; 1.5 bar system pressure.
- Not found: a primary-source flow-per-kW figure; case or junction temperature limits (not searched).

## Reduced-order models
| Model | Agreement | Checked against | Source |
|---|---|---|---|
| Manifold microchannel correlations, fit to 490 CFD cases | Friction within 0.7%; Nu mean 0.4%, worst 11–12% | CFD only | Yuruker, Mandel, Shooshtari, Ohadi, ITherm 2019 / Electronics Cooling 2022 |
| Porous-zone manifold model | Pressure drop within 5% | Experiment; thermal side not validated (*snippet*) | 2024 conference paper |
| Porous-medium microchannel layer | Not retrieved | 3D single-channel CFD (*snippet*) | Nonino & Savino 2025 |

## Sources read in full or in part
- Copeland, "Review of Low Profile Cold Plate Technology for High Density Servers", Electronics Cooling, 2005: https://www.electronics-cooling.com/2005/05/review-of-low-profile-cold-plate-technology-for-high-density-servers/
- Yuruker et al., "Heat Transfer and Pressure Drop Correlations for Manifold Microchannel Heatsinks", Electronics Cooling, 2022: https://www.electronics-cooling.com/2022/10/heat-transfer-and-pressure-drop-correlations-for-manifold-microchannel-heat-sinks/
- Robinson, IEEE Trans. CPMT 32(2), 2009 (abstract): https://edepositireland.ie/handle/2262/31928
- Lee & Vafai, Int. J. Heat Mass Transfer 42(9), 1999 (abstract): https://portal.mardi4nfdi.de/wiki/Item:Q1973216
- Fibermall GB200 cold plate listing (vendor): https://www.fibermall.com/blog/nvidia-gb200-liquid-cooling-plate.htm
- QCT Qoolrack white paper for GB200 NVL72 (vendor): https://blog.qct.io/wp-content/uploads/2025/04/QCT-Qoolrack-Stand-Alone_Advanced-Liquid-Cooling-for-NVIDIA-GB200-NVL72-Systems.pdf
- Kong et al., J. Electron. Packag. 2025 (snippet): https://atsla.cau.ac.kr/wp-content/uploads/2025/09/73-Kong-et-al-JEP2025.pdf
- OCP Cold Plate Test Specification V1.0.0 (unreadable PDF): https://www.opencompute.org/wp-content/uploads/2026/10/OCP-Cold-Plate-Test-Specification-V1.0.0.pdf

## Follow-up
Read the full text of Kong 2024/2025, Jung 2020 and Robinson 2009; obtain the OCP documents in readable form; find a primary
NVIDIA or partner coolant specification; search for accelerator case-temperature limits and for equal-pumping-power comparisons.
