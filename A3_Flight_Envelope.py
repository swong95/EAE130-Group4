# Flight Envelope

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Tuple

G = 32.174          # ft/s^2, gravitational acceleration
FT_PER_NMI = 6076.1 # ft per nautical mile

def atmosphere(h: float) -> Tuple[float, float]:
    """
    Standard atmosphere up to 65,000 ft (troposphere + lower stratosphere).

    Parameters
    ----------
    h : float
        Altitude in feet.

    Returns
    -------
    rho : float
        Air density [slug/ft^3]
    a   : float
        Speed of sound [ft/s]
    """
    # Sea-level standard values
    T0   = 518.67   # Rankine
    P0   = 2116.22  # lbf/ft^2
    rho0 = 0.002377 # slug/ft^3

    if h <= 36089:                          # Troposphere
        T   = T0 - 0.003566 * h            # Rankine
        rat = T / T0
        P   = P0 * rat ** 5.2561
        rho = rho0 * rat ** 4.2561
    else:                                   # Stratosphere (isothermal)
        T   = 389.97                        # Rankine (constant)
        P   = 472.68 * np.exp(-4.8063e-5 * (h - 36089))
        rho = P / (1716.49 * T)

    a = np.sqrt(1.4 * 1716.49 * T)         # ft/s
    return rho, a


@dataclass
class AircraftParams:
    """All fixed design parameters for the sizing analysis."""

    # --- Geometry ---
    S_ref: float = 600.0        # Wing reference area [ft^2]
    AR:    float = 2.08815          # Aspect ratio [-]  (F/A-18C ~= 3.5)
    e:     float = 0.80         # Oswald efficiency factor [-]
    tc:    float = 0.09         # Thickness-to-chord ratio at root [-]
    taper: float = 0.30         # Wing taper ratio lambda [-]
    sweep: float = 40.0         # Quarter-chord sweep [deg]

    # --- Aerodynamics ---
    CD0:   float = 0.01070       # Zero-lift drag coefficient [-]
    # k is derived from AR and e below

    # --- Propulsion (single engine) ---
    T_sl_max:  float = 43000  # Max sea-level static thrust per engine [lbf]
    
    n_engines: int   = 1
    ct_cruise: float = 0.866     # TSFC at cruise [1/hr]
    ct_loiter: float = 0.80     # TSFC at loiter [1/hr]
    ct_sl:     float = 0.40     # TSFC at sea-level idle (taxi) [1/hr]
    idle_frac: float = 0.05     # Idle thrust as fraction of max T/W

    # Thrust lapse with altitude: T = T_sl * (rho/rho0)^n
    thrust_lapse_exp: float = 0.7

    # Thrust scales with W0 to hold T/W constant (the engine grows with the aircraft).
    # T/W_target is derived from the PDR value: T_sl_max * n_engines / W_dg_PDR.
    # We store the PDR design weight as the reference point.
    W_dg_PDR:   float = 43200.0   # lbf  — PDR gross weight, used to fix T/W ratio
    TW_target:  float = field(init=False)  # derived in __post_init__

    # --- Weights (fixed) ---
    W_payload_a2a:    float = 5020   # Payload for air-to-air loadout [lbf]
    W_payload_strike: float = 6932   # Payload for ground strike loadout [lbf]

    # -----------------------------------------------------------------------
    # --- Component weight parameters (Raymer Ch. 15 + VSP geometry) ---
    # -----------------------------------------------------------------------

    # --- Raymer K-factors ---
    K_cb:  float = 1.0      # main landing gear crossbeam value
    K_d:   float = 2.75     # duct constant
    K_dw:  float = 0.768    # for delta wing
    K_dwf: float = 0.774    # for delta wing fuselage
    K_mc:  float = 1.45     # electronics mission completion
    K_rht: float = 1.047    # for rolling tail (2 surfaces move independently)
    K_tpg: float = 1.0      # tripod gear
    K_vg:  float = 1.62     # variable geometry air induction
    K_vs:  float = 1.0      # variable sweep wing
    K_vsh: float = 1.0      # variable sweep wing hydraulics

    # --- Lengths / areas (VSP geometry) ---
    L_d:   float = 5.809    # ft  duct length
    L_ec:  float = 10.503   # ft  engine front to cockpit
    L_sh:  float = 220/12   # ft  shroud length
    L_t:   float = 8.72490  # ft  tail length
    L_m:   float = 60.0     # in  main gear length
    L_n:   float = 60 + (2.704*12)  # in  nose gear length
    S_fw:  float = 400.0    # ft^2 firewall surface area
    S_vt:  float = 100.0    # ft^2 vertical tail area
    S_r:   float = 100.0    # ft^2 rudder area
    S_csw: float = 200.0    # ft^2 wing control surface area
    H_tH_v: float = 0.0     # no horizontal tail

    # --- Fuselage geometry ---
    L_fus: float = 40.28139 # ft  fuselage structural length
    D_fus: float = 5.4      # ft  fuselage depth
    W_fus: float = 13.232   # ft  fuselage structural width
    D_e:   float = 43/12    # ft  engine inlet diameter
    V_pr:  float = 120.0    # ft^3 pressurised section volume

    # --- Tail geometry ---
    AR_vt:    float = 1.00176
    tc_vt:    float = 0.100
    taper_vt: float = 0.25   # c_tip_vt / c_root_vt = 2.29229 / 9.16916
    sweep_vt: float = 45.0   # deg, 25%-chord sweep

    # --- Fuel / tanks ---
    V_i:  float = 120 * 7.48052  # gal  integral tank volume
    N_t:  float = 1.0             # number of fuel tanks
    SFC:  float = 0.886           # lb/lb/hr (used only in fuel system weight eq.)

    # --- Propulsion / engine ---
    W_en:  float = 6422.0   # lb  uninstalled engine weight (per engine)

    # --- Systems ---
    N_c:   float = 1.0      # number of crew
    N_ci:  float = 1.0      # crew config factor (1.0 = single pilot)
    N_s:   float = 1.0      # number of flight control systems
    N_u:   float = 10.0     # number of hydraulic utility functions
    N_nw:  float = 2.0      # number of nose wheels
    N_gear: float = 3.0     # total number of landing gear legs
    R_kva: float = 160.0    # kVA  electrical system rating
    W_uav: float = 2500.0   # lb   uninstalled avionics weight

    # --- Fudge factors (composite construction) ---
    wing_fudge:     float = 0.85
    tail_fudge:     float = 0.83
    fuselage_fudge: float = 0.90
    AIS_fudge:      float = 0.85

    # --- Nz / loads ---
    Nz: float = 10.5        # ultimate load factor (1.5 * 7g limit load)
    limloadf: float = 7.0   # limit load factor [g]

    # --- Mission profile ---
    # Climb
    h_cruise:     float = 40000.0   # Cruise altitude [ft]
    h_climb_step: float = 1000.0    # Altitude step size for climb integration [ft]

    # Cruise
    M_cruise:     float = 0.85      # Cruise Mach number [-]
    R_cruise:     float = 700*2    # Cruise range [nmi]
    n_cruise_seg: int   = 20        # Number of cruise segments

    # Combat / loiter
    E_loiter:     float = 0.333      # Loiter endurance [hr]

    # --- PDR reference weight — used to fix T/W through the sizing loop ---
    W_dg_PDR: float = 43200.0   # lbf, PDR gross weight estimate

    # --- Derived properties (computed in __post_init__) ---
    k:         float = field(init=False)
    TW_target: float = field(init=False)

    def __post_init__(self):
        self.k         = 1.0 / (np.pi * self.AR * self.e)
        # Lock in the T/W ratio from the PDR so the engine scales with the aircraft
        self.TW_target = (self.T_sl_max * self.n_engines) / self.W_dg_PDR

    def T_sl_for_W0(self, W0: float) -> float:
        """Sea-level static thrust scaled to maintain PDR T/W ratio [lbf]."""
        return self.TW_target * W0


# ---------------------------------------------------------------------------
# 3.  AERODYNAMIC HELPERS
# ---------------------------------------------------------------------------

def CL_from_weight(W: float, rho: float, V: float, p: AircraftParams) -> float:
    """Lift coefficient required to support weight W at speed V, density rho."""
    q = 0.5 * rho * V**2
    return W / (q * p.S_ref)


def CD_from_CL(CL: float, p: AircraftParams) -> float:
    """Quadratic drag polar."""
    return p.CD0 + p.k * CL**2


def LD(CL: float, p: AircraftParams) -> float:
    """Lift-to-drag ratio."""
    return CL / CD_from_CL(CL, p)


def LD_max(p: AircraftParams) -> float:
    """Maximum lift-to-drag ratio for the quadratic polar."""
    return 1.0 / (2.0 * np.sqrt(p.CD0 * p.k))


# ---------------------------------------------------------------------------
# 4.  THRUST LAPSE MODEL
# ---------------------------------------------------------------------------

def thrust_available(h: float, p: AircraftParams, W0: float = None) -> float:
    """
    Available max continuous thrust at altitude h [lbf].
    Uses simplified density-ratio lapse: T = T_sl * (rho/rho0)^n

    If W0 is supplied, sea-level thrust is scaled to maintain the PDR T/W
    ratio across sizing iterations (i.e. the engine grows with the aircraft).
    If W0 is None, uses the fixed T_sl_max * n_engines from params.
    """
    rho,  _ = atmosphere(h)
    rho0, _ = atmosphere(0.0)
    T_sl = p.T_sl_for_W0(W0) if W0 is not None else p.T_max_total
    return T_sl * (rho / rho0) ** p.thrust_lapse_exp

ceiling = 40000
alt = range(0, ceiling, 1000)

rho_arr = []
a_arr = []

for h in alt:
    rho, a = atmosphere(h)
    rho_arr.append(rho)
    a_arr.append(a)
    
rho_arr = np.array(rho_arr)
a_arr = np.array(a_arr)

w_strike = 52000
stall_arr = ((2*52000)/(rho_arr*AircraftParams.S_ref*2.03))**0.5



# print('stall_arr', stall_arr)
plt.figure(1)
plt.plot(stall_arr, alt)
plt.show()