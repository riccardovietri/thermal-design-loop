"""Typed, JSON-serialisable inputs and outputs shared by the tools."""
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Fluid:
    name: str
    temp_c: float
    rho: float  # kg/m3
    cp: float  # J/kg.K
    k: float  # W/m.K
    mu: float  # Pa.s


@dataclass(frozen=True)
class Requirements:
    heat_load_w: float
    footprint_width_m: float
    flow_length_m: float
    base_thickness_m: float
    solid_k: float
    inlet_temp_c: float
    dp_budget_pa: float
    coolant: str
    r_th_threshold_k_per_w: float | None = None  # must meet (SPEC Amendment 3)
    r_th_objective_k_per_w: float | None = None  # want
    marginal_return_cutoff: float = 0.25


@dataclass(frozen=True)
class Design:
    w_ch_m: float
    t_fin_m: float
    h_fin_m: float
    flow_lpm: float

    @property
    def design_id(self) -> str:
        return f"w{self.w_ch_m*1e6:04.0f}_t{self.t_fin_m*1e6:04.0f}_h{self.h_fin_m*1e6:04.0f}_q{self.flow_lpm*100:03.0f}"

    def n_channels(self, footprint_width_m: float) -> int:
        return int((footprint_width_m + 1e-12) // (self.w_ch_m + self.t_fin_m))


@dataclass(frozen=True)
class MeshInfo:
    path: str
    n_nodes: int
    n_elements: int
    n_fluid_elements: int
    min_quality_sicn: float
    mean_quality_sicn: float
    size_m: float


@dataclass(frozen=True)
class CfdResult:
    design_id: str
    n_channels: int
    reynolds: float
    dp_pa: float
    pump_power_w: float
    f_re_fanning: float
    nu_eff: float
    r_th_k_per_w: float
    r_th_max_k_per_w: float
    base_avg_temp_c: float
    base_max_temp_c: float
    coolant_rise_k: float
    energy_residual_rel: float
    n_dofs: int
    linear_residual_rel: float = 0.0
    base_flux_error_rel: float = 0.0


def to_dict(obj) -> dict:
    return asdict(obj)


@dataclass(frozen=True)
class DesignSpace:
    w_ch_m: tuple[float, ...]
    t_fin_m: tuple[float, ...]
    h_fin_m: tuple[float, ...]
    flow_lpm: tuple[float, ...]
    min_feature_m: float
    max_aspect_ratio: float


@dataclass(frozen=True)
class Screen1DResult:
    design: Design
    n_channels: int
    reynolds: float
    dp_pa: float
    r_th_k_per_w: float
    r_th_max_k_per_w: float
    r_base: float
    r_conv: float
    r_caloric: float
    fin_efficiency: float
    nusselt: float
    f_re_fanning: float
    feasible: bool
    violations: tuple[str, ...]


@dataclass(frozen=True)
class Screen1DEntranceResult(Screen1DResult):
    """Screen1DResult with developing-flow and inlet/outlet terms (`analyze_1d(..., entrance_effects=True)`).

    dp_pa is the channel-core pressure drop with hydrodynamically developing flow; r_th_k_per_w uses the simultaneously
    developing (combined-entry) heat transfer. The fully developed values and each contribution are kept separately.
    """
    dp_fully_developed_pa: float = 0.0
    dp_hydrodynamic_entrance_pa: float = 0.0  # extra core pressure drop of developing flow (f_app vs f)
    dp_inlet_outlet_pa: float = 0.0  # abrupt contraction + abrupt expansion, outside the channel core
    dp_with_inlet_outlet_pa: float = 0.0  # dp_pa + dp_inlet_outlet_pa
    r_th_fully_developed_k_per_w: float = 0.0
    r_th_thermal_entrance_k_per_w: float = 0.0  # thermal entrance only (Graetz: developed velocity)
    nusselt_fully_developed: float = 0.0
    nu_factor_thermal: float = 1.0  # effective h / developed h, thermal entrance only
    nu_factor_combined: float = 1.0  # effective h / developed h, simultaneously developing
    f_app_factor: float = 1.0  # apparent / developed friction factor over the channel length
    k_contraction: float = 0.0
    k_expansion: float = 0.0
    sigma: float = 0.0  # free-flow / frontal area, w / (w + t)
    x_star: float = 0.0  # L / (Dh Re Pr)
    l_plus_dh: float = 0.0  # L / (Dh Re)
