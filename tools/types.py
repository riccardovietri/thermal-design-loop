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


def to_dict(obj) -> dict:
    return asdict(obj)
