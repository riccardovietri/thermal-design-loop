"""OpenFOAM case for the streamwise-periodic conjugate unit cell, run in Docker, and parsers for its output.

Periodic formulation (Patankar, Liu & Sparrow 1977), the same as Level 1's fully developed one:
  flow:   cyclic in x, the mean velocity held at the Level 1 value by `meanVelocityForce`; the pressure gradient it needs is
          the fully developed dp/dx.
  energy: T = gamma x + theta(x, y, z) with theta periodic; solving for theta adds the sink -rho cp gamma u_x in the fluid
          (nothing in the solid: the Laplacian of gamma x is zero). The sink is written as -Q_cell u_x / sum(u_x V), which is
          the same field at convergence and removes all heat that enters the base, so the theta level cannot drift.
"""
import os
import re
import subprocess
from pathlib import Path

from .geometry import UnitCell, _header, blockmesh_dict, grid

# Pinned by digest: OpenCFD OpenFOAM v2406, the "default" image (includes wmake, needed for the coded sink).
IMAGE = "opencfd/openfoam-default:2406@sha256:dd5aa20630a55722663bf83ba0cb74870cba130081303e32e3865007fa2aa35a"
T_REF_K = 313.15  # theta level: the solve starts at the inlet temperature; only differences are used
MOLAR_MASS = 18.0  # rhoConst / hConst: irrelevant to the solution, required by the thermophysical library
CU_RHO, CU_CP = 8960.0, 385.0  # steady state: neither enters the solution


def _w(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _field(obj: str, cls: str, dims: str, internal: str, patches: dict[str, str]) -> str:
    body = "".join(f"    {p}\n    {{\n{spec}    }}\n" for p, spec in patches.items())
    return (_header(obj, cls) + f"dimensions {dims};\n\ninternalField uniform {internal};\n\nboundaryField\n{{\n{body}}}\n")


def _p(**kv) -> str:
    return "".join(f"        {k} {v};\n" for k, v in kv.items())


def write_case(case: Path, cell: UnitCell, factor: int, n_iter: int) -> dict:
    case = Path(case)
    g = grid(cell, factor)
    q_cell = cell.cell_heat_w(g["lx_m"])
    for region in ("fluid", "solid"):
        _w(case / "system" / region / "blockMeshDict", blockmesh_dict(cell, region, factor))
    _w(case / "constant/regionProperties", _header("regionProperties") + "regions (fluid (fluid) solid (solid));\n")
    # --- fluid ---
    pr = cell.cp * cell.mu / cell.k_fluid
    _w(case / "constant/fluid/thermophysicalProperties", _header("thermophysicalProperties") + f"""
thermoType {{ type heRhoThermo; mixture pureMixture; transport const; thermo hConst; equationOfState rhoConst;
             specie specie; energy sensibleEnthalpy; }}
mixture
{{
    specie {{ molWeight {MOLAR_MASS}; }}
    equationOfState {{ rho {cell.rho!r}; }}
    thermodynamics {{ Cp {cell.cp!r}; Hf 0; }}
    transport {{ mu {cell.mu!r}; Pr {pr!r}; }}
}}
""")
    _w(case / "constant/fluid/turbulenceProperties", _header("turbulenceProperties") + "simulationType laminar;\n")
    _w(case / "constant/fluid/radiationProperties", _header("radiationProperties") + "radiation off;\nradiationModel none;\n")
    _w(case / "constant/g", _header("g", "uniformDimensionedVectorField") + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    _w(case / "constant/fluid/fvOptions", _header("fvOptions") + f"""
meanFlow
{{
    type            meanVelocityForce;
    selectionMode   all;
    fields          (U);
    Ubar            ({cell.u_bar!r} 0 0);
}}

periodicSink
{{
    type            scalarCodedSource;
    selectionMode   all;
    fields          (h);
    name            periodicSink;
    codeInclude
    #{{
    #}};
    codeCorrect
    #{{
    #}};
    codeConstrain
    #{{
    #}};
    codeAddSup
    #{{
        FatalErrorInFunction << "compressible solver expected" << exit(FatalError);
    #}};
    codeAddSupRho
    #{{
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
        const scalarField& V = mesh_.V();
        const scalarField ux(U.primitiveField().component(0));
        const scalar sumUV = gSum(ux*V);
        const scalar Qcell = {q_cell!r};
        eqn.source() += Qcell*ux*V/sumUV;
        Info<< "L3MON bulkT " << gSum(ux*T.primitiveField()*V)/sumUV << " flow " << sumUV/{g['lx_m']!r}
            << " sink " << Qcell << endl;
    #}};
}}
""")
    wall = {"U": "        type noSlip;\n", "T": "        type zeroGradient;\n",
            "p_rgh": "        type fixedFluxPressure;\n        value uniform 1e5;\n", "p": "        type calculated;\n        value uniform 1e5;\n"}
    coupled = ("        type compressible::turbulentTemperatureCoupledBaffleMixed;\n"
               "        Tnbr T;\n        kappaMethod {km};\n        value uniform " + f"{T_REF_K};\n")
    cyc = "        type cyclic;\n"
    sym = "        type symmetryPlane;\n"
    fpat = lambda f, ifc: {"front": cyc, "back": cyc, "symFluid": sym, "cover": wall[f], "fluid_to_solid": ifc}
    _w(case / "0/fluid/U", _field("U", "volVectorField", "[0 1 -1 0 0 0 0]", f"({cell.u_bar!r} 0 0)", fpat("U", wall["U"])))
    _w(case / "0/fluid/T", _field("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{T_REF_K}",
                                  fpat("T", coupled.format(km="fluidThermo"))))
    _w(case / "0/fluid/p_rgh", _field("p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5", fpat("p_rgh", wall["p_rgh"])))
    _w(case / "0/fluid/p", _field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5", fpat("p", wall["p"])))
    # --- solid ---
    _w(case / "constant/solid/thermophysicalProperties", _header("thermophysicalProperties") + f"""
thermoType {{ type heSolidThermo; mixture pureMixture; transport constIso; thermo hConst; equationOfState rhoConst;
             specie specie; energy sensibleEnthalpy; }}
mixture
{{
    specie {{ molWeight 63.5; }}
    transport {{ kappa {cell.k_solid!r}; }}
    thermodynamics {{ Hf 0; Cp {CU_CP}; }}
    equationOfState {{ rho {CU_RHO}; }}
}}
""")
    _w(case / "constant/solid/radiationProperties", _header("radiationProperties") + "radiation off;\nradiationModel none;\n")
    flux = ("        type externalWallHeatFluxTemperature;\n        mode flux;\n"
            f"        q uniform {cell.heated_flux_w_m2!r};\n        kappaMethod solidThermo;\n        value uniform {T_REF_K};\n")
    spat = {"front": cyc, "back": cyc, "symFin": sym, "symBase": sym, "finTop": "        type zeroGradient;\n", "bottom": flux,
            "solid_to_fluid": coupled.format(km="solidThermo")}
    _w(case / "0/solid/T", _field("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{T_REF_K}", spat))
    _w(case / "0/solid/p", _field("p", "volScalarField", "[1 -1 -2 0 0 0 0]", "1e5",
                                  {k: (v if v in (cyc, sym) else "        type calculated;\n        value uniform 1e5;\n")
                                   for k, v in spat.items()}))
    # --- numerics ---
    _w(case / "system/fluid/fvSchemes", _header("fvSchemes") + """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div(phi,h) bounded Gauss linearUpwind grad(h);
    div(phi,K) bounded Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    _w(case / "system/fluid/fvSolution", _header("fvSolution") + f"""
solvers
{{
    rho {{ solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }}
    p_rgh {{ solver GAMG; smoother GaussSeidel; tolerance 1e-12; relTol 0.01; maxIter 20; }}
    "(U|h)" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-14; relTol 0.01; }}
}}
SIMPLE
{{
    momentumPredictor yes;
    nNonOrthogonalCorrectors 0;
    pRefCell 0;
    pRefValue 1e5;
}}
relaxationFactors
{{
    fields {{ rho 1; p_rgh 0.7; }}
    equations {{ U 0.7; h 0.9; }}
}}
""")
    _w(case / "system/solid/fvSchemes", _header("fvSchemes") + """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    _w(case / "system/solid/fvSolution", _header("fvSolution") + f"""
solvers {{ h {{ solver PCG; preconditioner DIC; tolerance 1e-14; relTol 0.01; }} }}
SIMPLE {{ nNonOrthogonalCorrectors 0; }}
relaxationFactors {{ equations {{ h 1; }} }}
""")
    _w(case / "system/fvSchemes", _header("fvSchemes"))
    _w(case / "system/fvSolution", _header("fvSolution"))
    _w(case / "system/controlDict", _header("controlDict") + f"""
application     chtMultiRegionSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {n_iter};
deltaT          1;
writeControl    timeStep;
writeInterval   {n_iter};
purgeWrite      1;
writeFormat     ascii;
writePrecision  12;
timePrecision   8;
runTimeModifiable false;

functions
{{
    baseT
    {{
        type surfaceFieldValue; libs (fieldFunctionObjects); region solid; regionType patch; name bottom;
        operation areaAverage; fields (T); writeFields false; writeControl timeStep; writeInterval 1; log false;
    }}
    residualsFluid
    {{
        type solverInfo; libs (utilityFunctionObjects); region fluid; fields (U p_rgh h); writeResidualFields no;
    }}
    residualsSolid
    {{
        type solverInfo; libs (utilityFunctionObjects); region solid; fields (h); writeResidualFields no;
    }}
    heatFluid
    {{
        type wallHeatFlux; libs (fieldFunctionObjects); region fluid; patches (fluid_to_solid);
        writeControl writeTime; writeFields false;
    }}
    heatSolid
    {{
        type wallHeatFlux; libs (fieldFunctionObjects); region solid; patches (bottom solid_to_fluid);
        writeControl writeTime; writeFields false;
    }}
}}
""")
    return {"grid": g, "q_cell_w": q_cell}


ALLRUN = """set -e
cd /case
blockMesh -region fluid > log.blockMesh.fluid 2>&1
blockMesh -region solid > log.blockMesh.solid 2>&1
checkMesh -region fluid > log.checkMesh.fluid 2>&1
checkMesh -region solid > log.checkMesh.solid 2>&1
chtMultiRegionSimpleFoam > log.solver 2>&1 || true
postProcess -region fluid -func writeCellCentres -latestTime > log.centres.fluid 2>&1
postProcess -region solid -func writeCellCentres -latestTime > log.centres.solid 2>&1
postProcess -region fluid -func writeCellVolumes -latestTime > log.volumes.fluid 2>&1
cat /sys/fs/cgroup/memory.peak > memory.peak
"""


def run_docker(case: Path) -> None:
    case = Path(case).resolve()
    _w(case / "Allrun", ALLRUN)
    subprocess.run(["docker", "run", "--rm", "--cpus", "1", "--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp", "-v", f"{case}:/case", "-w", "/case", IMAGE, "-c", "bash /case/Allrun"],
                   check=True)


# --- parsers ---

def parse_solver_log(text: str) -> dict:
    """Per-iteration monitored values from log.solver: bulk T (coded sink), mean-velocity-force pressure gradient."""
    bulk = [float(m.group(1)) for m in re.finditer(r"L3MON bulkT (\S+) flow (\S+)", text)]
    flow = [float(m.group(2)) for m in re.finditer(r"L3MON bulkT (\S+) flow (\S+)", text)]
    blocks = re.split(r"^Time = \S+\s*$", text, flags=re.M)[1:]
    # dp/dx held by meanVelocityForce: the last value printed in each iteration (after the pressure corrector).
    # For the compressible solver it is in Pa/m (the force is per unit volume), checked against the exact duct f·Re.
    grad = [float(re.findall(r"[Pp]ressure gradient = (\S+)", b)[-1]) for b in blocks if "ressure gradient" in b]
    times = [int(float(m.group(1))) for m in re.finditer(r"^Time = (\S+)", text, re.M)]
    return {"iterations": times, "bulk_t": bulk, "flow_m3_s": flow, "dpdx_pa_per_m": grad,
            "converged": "reached convergence" in text or "SIMPLE solution converged" in text,
            "ended": "End" in text.split("\n")[-3:] or bool(re.search(r"^End\s*$", text, re.M))}


def parse_dat(path: Path) -> list[list[float]]:
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        vals = []
        for tok in line.replace("(", " ").replace(")", " ").split():
            try:
                vals.append(float(tok))
            except ValueError:
                vals.append(float("nan"))
        rows.append(vals)
    return rows


def parse_check_mesh(text: str) -> dict:
    def num(pattern):
        m = re.search(pattern, text)
        return float(m.group(1)) if m else None
    return {"cells": num(r"cells:\s+(\d+)"), "max_non_orthogonality": num(r"Mesh non-orthogonality Max: (\S+)"),
            "max_skewness": num(r"Max skewness = (\S+)"), "max_aspect_ratio": num(r"Max aspect ratio = (\S+)"),
            "mesh_ok": "Mesh OK." in text}
