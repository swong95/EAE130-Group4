"""
Refined Aircraft Sizing - F/A-18 Replacement
=============================================
Implements the iterative weight estimation loop described in B-02-SizingRefinement,
using segmented climb (energy-height method) and segmented exponential Breguet cruise.

All units are imperial:
    Weight / Force : lbf
    Length / Alt   : ft
    Speed          : ft/s  (converted to knots or nmi where noted)
    Distance       : nmi   (nautical miles)
    Time           : hr
    Density        : slug/ft^3
    TSFC (ct)      : 1/hr
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Tuple

# ---------------------------------------------------------------------------
# 1.  CONSTANTS & ATMOSPHERE
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# 2.  AIRCRAFT DESIGN PARAMETERS  (edit these to match your PDR values)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# 5.  BEST CLIMB SPEED FOR A JET  (slide 48)
# ---------------------------------------------------------------------------

def best_climb_speed(W: float, h: float, p: AircraftParams,
                     W0: float = None) -> float:
    """
    Speed [ft/s] that maximises rate of climb for a jet at altitude h,
    weight W, assuming thrust = T_available (max continuous).

    W0 is the current sizing-loop gross weight, used to scale thrust
    so T/W stays constant. Pass None to use fixed T_sl_max from params.
    """
    rho, _ = atmosphere(h)
    T      = thrust_available(h, p, W0=W0)
    TW     = T / W
    WS     = W / p.S_ref
    A      = TW * WS / rho
    B      = 12.0 * p.CD0 * p.k * (WS / rho) ** 2
    V2     = A + np.sqrt(A**2 + B)
    return np.sqrt(max(V2, 1.0))


# ---------------------------------------------------------------------------
# 6.  FUEL FRACTIONS  — TAXI / TAKEOFF
# ---------------------------------------------------------------------------

def fuel_fraction_taxi(p: AircraftParams) -> float:
    """
    Taxi + startup: 15 min at idle thrust, sea-level TSFC.
    W_{i+1}/W_i = 1 - t * ct * (T_idle / W)
    We fold the T/W ratio into the calculation at sea level.
    """
    t_hr    = 15.0 / 60.0                  # 15 min -> hours
    T_idle  = p.T_max_total * p.idle_frac  # lbf at idle
    # ct_sl in 1/hr, T_idle/W gives dimensionless thrust fraction
    # We return the weight ratio; W cancels because we track ratios
    # dW/dt = -ct * T  =>  W_{f}/W_{i} = 1 - ct*t*(T/W) * 1
    # But T/W depends on W, which we don't know yet.
    # Approximation: use a representative T/W at idle ~ (T_idle / W0_guess)
    # This is small enough that the approximation is fine.
    # We return the *fractional fuel burn* per unit W; caller scales.
    fuel_frac_idle = t_hr * p.ct_sl * p.idle_frac   # unitless approximation
    return max(0.0, fuel_frac_idle)                   # fuel burned / W_start


def fuel_fraction_takeoff(p: AircraftParams) -> float:
    """
    Takeoff: 1 min at max thrust, sea-level TSFC.
    Fuel fraction burned = ct * t * (T/W); T/W ~ 1 at max thrust.
    """
    t_hr = 1.0 / 60.0
    return p.ct_sl * t_hr * 1.0     # T/W = 1 at max thrust (approx)


# ---------------------------------------------------------------------------
# 7.  CLIMB SEGMENT  (energy-height method, slide 47)
# ---------------------------------------------------------------------------

def simulate_climb(W_start: float, p: AircraftParams,
                   W0: float = None) -> dict:
    """
    Simulate climb from sea level to h_cruise using the energy-height method.

    W0 is the current sizing-loop gross weight passed in so that thrust
    scales correctly with the aircraft size each iteration.
    """
    h      = 0.0
    W      = W_start
    h_end  = p.h_cruise
    dh     = p.h_climb_step

    hist = {
        'h': [], 'V_kts': [], 'W': [], 'RC': [], 'gamma_deg': [],
        'fuel_cum': [], 'x_nmi': [], 'LD': []
    }

    fuel_total = 0.0
    x_total    = 0.0

    while h < h_end:
        dh_step = min(dh, h_end - h)
        h_mid   = h + 0.5 * dh_step

        rho_mid, _ = atmosphere(h_mid)
        T_mid      = thrust_available(h_mid, p, W0=W0)

        V   = best_climb_speed(W, h_mid, p, W0=W0)
        CL  = CL_from_weight(W, rho_mid, V, p)
        CD  = CD_from_CL(CL, p)
        D   = 0.5 * rho_mid * V**2 * p.S_ref * CD
        L_D = CL / CD

        Ps    = V * (T_mid - D) / W
        Ps    = max(Ps, 1.0)
        sin_g = Ps / V
        gamma = np.arcsin(np.clip(sin_g, -1, 1))

        delta_he = dh_step          # ≈ he_top - he_bot for small steps
        dt       = delta_he / Ps    # seconds

        ct_climb_s = p.ct_cruise / 3600.0
        delta_W    = ct_climb_s * T_mid * dt
        delta_x    = V * np.cos(gamma) * dt

        fuel_total += delta_W
        x_total    += delta_x
        W          -= delta_W
        h          += dh_step

        hist['h'].append(h_mid)
        hist['V_kts'].append(V / 1.6878)
        hist['W'].append(W)
        hist['RC'].append(Ps * 60)
        hist['gamma_deg'].append(np.degrees(gamma))
        hist['fuel_cum'].append(fuel_total)
        hist['x_nmi'].append(x_total / FT_PER_NMI)
        hist['LD'].append(L_D)

    return {
        'W_end':       W,
        'fuel_burned': fuel_total,
        'x_climb':     x_total / FT_PER_NMI,
        'history':     {k: np.array(v) for k, v in hist.items()}
    }


# ---------------------------------------------------------------------------
# 8.  CRUISE SEGMENT  (segmented exponential Breguet, slide 34)
# ---------------------------------------------------------------------------

def simulate_cruise(W_start: float, R_cruise_nmi: float, p: AircraftParams) -> dict:
    """
    Simulate cruise using segmented exponential Breguet.

    At each range segment ΔR:
        CL_i  = 2*W_i / (rho * V^2 * S)
        CD_i  = CD0 + k*CL_i^2
        LD_i  = CL_i / CD_i
        W_{i+1} = W_i * exp( -ΔR * ct / (V * LD_i) )

    All range in ft internally, converted from nmi at start.

    Returns
    -------
    dict with keys:
        W_end        : weight at end of cruise [lbf]
        fuel_burned  : total fuel burned [lbf]
        history      : dict of arrays for plotting
    """
    rho, a = atmosphere(p.h_cruise)
    V      = p.M_cruise * a                     # True airspeed [ft/s]

    R_ft   = R_cruise_nmi * FT_PER_NMI         # Total cruise range [ft]
    dR     = R_ft / p.n_cruise_seg             # Segment length [ft]

    W      = W_start
    ct_s   = p.ct_cruise                        # 1/hr (Breguet uses V in ft/s, R in ft,
                                                #  so we need ct in 1/s)
    ct_per_s = ct_s / 3600.0

    hist = {
        'R_nmi': [], 'W': [], 'CL': [], 'CD': [], 'LD': [],
        'fuel_cum': [], 'T_req': []
    }

    fuel_total = 0.0
    R_covered  = 0.0

    for _ in range(p.n_cruise_seg):
        CL   = CL_from_weight(W, rho, V, p)
        CD   = CD_from_CL(CL, p)
        L_D  = CL / CD
        T_req = W / L_D                          # Required thrust [lbf]

        # Exponential Breguet step
        W_new = W * np.exp(-dR * ct_per_s / (V * L_D))
        delta_W = W - W_new

        fuel_total += delta_W
        R_covered  += dR
        W           = W_new

        hist['R_nmi'].append(R_covered / FT_PER_NMI)
        hist['W'].append(W)
        hist['CL'].append(CL)
        hist['CD'].append(CD)
        hist['LD'].append(L_D)
        hist['fuel_cum'].append(fuel_total)
        hist['T_req'].append(T_req)

    return {
        'W_end':       W,
        'fuel_burned': fuel_total,
        'history':     {k: np.array(v) for k, v in hist.items()}
    }


# ---------------------------------------------------------------------------
# 9.  LOITER SEGMENT  (slide 51, Breguet endurance for jets)
# ---------------------------------------------------------------------------

def simulate_loiter(W_start: float, p: AircraftParams) -> dict:
    """
    Loiter at (L/D)_max for E_loiter hours.

    For jets:  W_{i+1}/W_i = exp( -E * ct / (L/D) )

    Fly at (L/D)_max by choosing CL = sqrt(CD0/k).
    """
    CL_opt = np.sqrt(p.CD0 / p.k)
    CD_opt = CD_from_CL(CL_opt, p)
    LD_opt = CL_opt / CD_opt

    W_end = W_start * np.exp(-p.E_loiter * p.ct_loiter / LD_opt)
    fuel  = W_start - W_end

    return {
        'W_end':       W_end,
        'fuel_burned': fuel,
        'LD_loiter':   LD_opt,
        'CL_loiter':   CL_opt
    }


# ---------------------------------------------------------------------------
# 10.  DESCENT + LANDING  (historical fractions, slide 52)
# ---------------------------------------------------------------------------

FRAC_DESCENT = 0.990    # W_{i+1}/W_i historical fraction
FRAC_LANDING = 0.995


def fuel_descent(W_start: float) -> Tuple[float, float]:
    """Returns (W_end, fuel_burned)."""
    W_end = W_start * FRAC_DESCENT
    return W_end, W_start - W_end


def fuel_landing(W_start: float) -> Tuple[float, float]:
    """Returns (W_end, fuel_burned)."""
    W_end = W_start * FRAC_LANDING
    return W_end, W_start - W_end


# ---------------------------------------------------------------------------
# 11.  EMPTY WEIGHT BUILDUP
# ---------------------------------------------------------------------------

def empty_weight(W0: float, p: AircraftParams) -> dict:
    """
    Full Raymer Ch. 15 component weight buildup.

    W0 is the current iteration's gross takeoff weight guess [lbf].
    All equations that previously used the fixed W_dg now use W0 so
    that the sizing loop is internally consistent.

    Returns a dict containing every component weight AND the total
    empty weight under the key 'W_empty', so the caller can inspect
    individual components or print a breakdown.
    """
    # Convenience aliases that match the Component_Weights_and_CG.py naming
    W_dg  = W0
    N_z   = p.Nz
    N_en  = p.n_engines
    # Scale thrust with W0 so component equations are consistent with sizing loop
    T     = p.T_sl_for_W0(W0)       # total sea-level thrust at current W0
    T_e   = T / N_en                 # per-engine thrust
    W_l   = W_dg * 0.78           # design landing gross weight
    N_l   = p.N_gear * 1.5        # ultimate landing load factor
    L_a   = p.L_ec                 # electrical routing distance

    # Fuel volume: keep proportional to W0 so the loop is self-consistent.
    # Base volume is p.V_i (set from PDR); we scale it if you later want
    # to make it a function of fuel weight — for now use the fixed value.
    V_i   = p.V_i
    V_t   = V_i
    V_p   = V_i

    sweep_rad    = np.radians(p.sweep)
    sweep_vt_rad = np.radians(p.sweep_vt)
    S_cs  = p.S_r + p.S_csw

    # ------------------------------------------------------------------
    # Structures group
    # ------------------------------------------------------------------
    W_wing = (0.0103 * p.K_dw * p.K_vs
              * (W_dg * N_z) ** 0.5
              * p.S_ref ** 0.622
              * p.AR ** 0.785
              * p.tc ** (-0.4)
              * (1 + p.taper) ** 0.05
              * np.cos(sweep_rad) ** (-1.0)
              * p.S_csw ** 0.04) * p.wing_fudge

    W_vt = (0.452 * p.K_rht
            * (1 + p.H_tH_v) ** 0.5
            * (W_dg * N_z) ** 0.488
            * p.S_vt ** 0.718
            * p.M_cruise ** 0.341
            * p.L_t ** (-1.0)
            * (1 + p.S_r / p.S_vt) ** 0.348
            * p.AR_vt ** 0.223
            * (1 + p.taper_vt) ** 0.25
            * np.cos(sweep_vt_rad) ** (-0.323)) * p.tail_fudge

    W_fuselage = (0.499 * p.K_dwf
                  * W_dg ** 0.35
                  * N_z ** 0.25
                  * p.L_fus ** 0.5
                  * p.D_fus ** 0.849
                  * p.W_fus ** 0.685) * p.fuselage_fudge

    W_mlg = (p.K_cb * p.K_tpg
             * (W_l * N_l) ** 0.25
             * p.L_m ** 0.973)

    W_nlg = ((W_l * N_l) ** 0.290
             * p.L_n ** 0.5
             * p.N_nw ** 0.525)

    W_ensect  = 0.01 * (p.W_en ** 0.717) * N_en * N_z
    W_handgear = 3.2e-4 * W_dg
    W_arrgear  = 0.008 * W_dg
    W_catgear  = 0.003 * W_dg
    W_canopy   = 200.0   # lb — fixed estimate based on F-22

    # ------------------------------------------------------------------
    # Propulsion group
    # ------------------------------------------------------------------
    W_emounts  = 0.013 * (N_en ** 0.795) * (T ** 0.579) * N_z
    W_firewall = 1.13 * p.S_fw
    W_encool   = 4.55 * p.D_e * p.L_sh * N_en
    W_oilcool  = 37.82 * (N_en ** 1.023)
    W_encontrols = 10.5 * (N_en ** 1.008) * (p.L_ec ** 0.222)
    W_starter  = 0.025 * (T_e ** 0.760) * (N_en ** 0.72)
    W_tailpipe = 3.5 * p.D_e * p.L_sh * N_en

    # L_s = single duct length = L_d for this design (ratio = 1.0)
    L_s_over_L_d = 1.0
    W_ais = (13.29 * p.K_vg
             * p.L_d ** 0.643
             * p.K_d ** 0.182
             * N_en ** 1.498
             * L_s_over_L_d ** (-0.373)
             * p.D_e) * p.AIS_fudge

    W_fuelsystanks = (7.45
                      * V_t ** 0.47
                      * (1 + V_i / V_t) ** (-0.095)
                      * (1 + V_p / V_t)
                      * p.N_t ** 0.066
                      * N_en ** 0.052
                      * (T * p.SFC / 1000) ** 0.249)

    # ------------------------------------------------------------------
    # Equipment / systems group
    # ------------------------------------------------------------------
    W_flightcont = (36.28
                    * p.M_cruise ** 0.003
                    * S_cs ** 0.489
                    * p.N_s ** 0.484
                    * p.N_c ** 0.127)

    W_instr = (8.0
               + 36.37 * N_en ** 0.676 * p.N_t ** 0.237
               + 26.4 * (1 + p.N_ci) ** 1.356)

    W_hydraulics = 37.23 * p.K_vsh * p.N_u ** 0.664

    W_elec = (172.2 * p.K_mc
              * p.R_kva ** 0.152
              * p.N_c ** 0.10
              * L_a ** 0.10
              * (p.n_engines ** 0.091))   # N_gen = N_en

    W_avionics = 2.117 * p.W_uav ** 0.933
    W_furn     = 217.6 * p.N_c
    W_acai     = 201.6 * ((p.W_uav + 200 * p.N_c) / 1000) ** 0.735

    # ------------------------------------------------------------------
    # Group totals
    # ------------------------------------------------------------------
    W_structures = (W_wing + W_vt + W_fuselage + W_ais + W_canopy
                    + W_mlg + W_nlg + W_ensect + W_handgear
                    + W_arrgear + W_catgear)

    W_propulsion = (p.W_en + W_fuelsystanks + W_emounts + W_encontrols
                    + W_encool + W_oilcool + W_tailpipe + W_starter)

    W_equipment  = (W_avionics + W_firewall + W_flightcont + W_instr
                    + W_hydraulics + W_elec + W_acai + W_furn)

    W_empty_total = W_structures + W_propulsion + W_equipment

    return {
        # --- Group totals ---
        'W_empty':       W_empty_total,
        'W_structures':  W_structures,
        'W_propulsion':  W_propulsion,
        'W_equipment':   W_equipment,
        # --- Individual components (for reporting) ---
        'W_wing':        W_wing,
        'W_vt':          W_vt,
        'W_fuselage':    W_fuselage,
        'W_mlg':         W_mlg,
        'W_nlg':         W_nlg,
        'W_ensect':      W_ensect,
        'W_handgear':    W_handgear,
        'W_arrgear':     W_arrgear,
        'W_catgear':     W_catgear,
        'W_canopy':      W_canopy,
        'W_ais':         W_ais,
        'W_en':          p.W_en,
        'W_fuelsystanks':W_fuelsystanks,
        'W_emounts':     W_emounts,
        'W_firewall':    W_firewall,
        'W_encool':      W_encool,
        'W_oilcool':     W_oilcool,
        'W_tailpipe':    W_tailpipe,
        'W_starter':     W_starter,
        'W_encontrols':  W_encontrols,
        'W_flightcont':  W_flightcont,
        'W_instr':       W_instr,
        'W_hydraulics':  W_hydraulics,
        'W_elec':        W_elec,
        'W_avionics':    W_avionics,
        'W_furn':        W_furn,
        'W_acai':        W_acai,
    }


# ---------------------------------------------------------------------------
# 12.  MAIN SIZING LOOP
# ---------------------------------------------------------------------------

def run_sizing(p: AircraftParams,
               W_payload: float = None,
               W0_guess: float = 45000.0,
               tol: float = 1.0,
               max_iter: int = 50,
               verbose: bool = True) -> dict:
    """
    Iterative weight estimation loop.

    Parameters
    ----------
    p         : AircraftParams
    W_payload : Payload weight to use [lbf]. If None, defaults to p.W_payload_a2a.
    W0_guess  : Initial gross takeoff weight guess [lbf]
    tol       : Convergence tolerance [lbf]
    max_iter  : Maximum iterations
    verbose   : Print iteration history

    Returns
    -------
    dict with converged weights, fuel breakdown, and simulation histories
    """
    if W_payload is None:
        W_payload = p.W_payload_a2a

    W0 = W0_guess
    converged = False

    # These will be filled on the last (converged) iteration
    last_climb  = None
    last_cruise = None
    last_loiter = None

    if verbose:
        print(f"\n{'='*60}")
        print(f"  Refined Sizing Loop  (tol = {tol:.1f} lbf)")
        print(f"{'='*60}")
        print(f"  {'Iter':>4}  {'W0_guess':>10}  {'W0_calc':>10}  {'Delta':>10}")
        print(f"  {'-'*44}")

    for iteration in range(1, max_iter + 1):

        # --- Phase 1: Taxi + Startup ---
        # Use scaled thrust (T/W constant) for fuel fractions
        T_sl_current = p.T_sl_for_W0(W0)
        ff_taxi = (15.0/60.0) * p.ct_sl * p.idle_frac          # idle frac of T/W
        ff_to   = (1.0 /60.0) * p.ct_sl * (T_sl_current / W0)  # full thrust 1 min

        W_after_taxi = W0 * (1 - ff_taxi)
        W_after_to   = W_after_taxi * (1 - ff_to)

        fuel_taxi = W0 - W_after_taxi
        fuel_to   = W_after_taxi - W_after_to

        # --- Phase 2: Climb ---
        climb = simulate_climb(W_after_to, p, W0=W0)
        W_after_climb = climb['W_end']
        fuel_climb    = climb['fuel_burned']
        x_climb       = climb['x_climb']        # nmi

        # --- Phase 3: Cruise ---
        R_remaining = max(0.0, p.R_cruise - x_climb)
        cruise = simulate_cruise(W_after_climb, R_remaining, p)
        W_after_cruise = cruise['W_end']
        fuel_cruise    = cruise['fuel_burned']

        # --- Phase 4: Combat Loiter ---
        loiter = simulate_loiter(W_after_cruise, p)
        W_after_loiter = loiter['W_end']
        fuel_loiter    = loiter['fuel_burned']

        # --- Phase 5: Descent ---
        W_after_desc, fuel_desc = fuel_descent(W_after_loiter)

        # --- Phase 6: Landing ---
        W_land, fuel_land = fuel_landing(W_after_desc)

        # --- Total fuel ---
        W_fuel_total = (fuel_taxi + fuel_to + fuel_climb
                        + fuel_cruise + fuel_loiter
                        + fuel_desc + fuel_land)

        # Add 5% reserve fuel
        W_fuel_total *= 1.05

        # --- Empty weight (full component buildup) ---
        ew      = empty_weight(W0, p)
        W_empty = ew['W_empty']

        # --- New W0 estimate ---
        W0_new = W_empty + W_fuel_total + W_payload

        delta = abs(W0_new - W0)

        if verbose:
            print(f"  {iteration:>4}  {W0:>10.1f}  {W0_new:>10.1f}  {delta:>10.2f}")

        # Save last iteration results
        last_climb  = climb
        last_cruise = cruise
        last_loiter = loiter

        if delta < tol:
            converged = True
            W0 = W0_new
            break

        W0 = W0_new

    if verbose:
        print(f"  {'-'*44}")
        if converged:
            print(f"  Converged in {iteration} iterations.")
        else:
            print(f"  WARNING: Did not converge in {max_iter} iterations.")
        print(f"{'='*60}\n")

    # Build fuel breakdown summary
    fuel_breakdown = {
        'Taxi & Startup': fuel_taxi,
        'Takeoff':        fuel_to,
        'Climb':          fuel_climb,
        'Cruise':         fuel_cruise,
        'Loiter':         fuel_loiter,
        'Descent':        fuel_desc,
        'Landing':        fuel_land,
    }
    fuel_reserve = W_fuel_total - sum(fuel_breakdown.values())
    fuel_breakdown['Reserve (5%)'] = fuel_reserve

    return {
        'converged':         converged,
        'W0':                W0,
        'W_empty':           W_empty,
        'W_fuel':            W_fuel_total,
        'W_payload':         W_payload,
        'fuel_breakdown':    fuel_breakdown,
        'component_weights': ew,
        'x_climb_nmi':       x_climb,
        'R_cruise_actual':   p.R_cruise - x_climb,
        'LD_max':            LD_max(p),
        'loiter_LD':         last_loiter['LD_loiter'],
        'climb_history':     last_climb['history'],
        'cruise_history':    last_cruise['history'],
    }


# ---------------------------------------------------------------------------
# 13.  RESULTS REPORTING
# ---------------------------------------------------------------------------

def print_results(res: dict, p: AircraftParams):
    """Pretty-print the sizing results including component weight breakdown."""
    cw = res['component_weights']
    print(f"\n{'='*60}")
    print(f"  SIZING RESULTS SUMMARY")
    print(f"{'='*60}")
    print(f"  Gross Takeoff Weight  W0   : {res['W0']:>10.1f} lbf")
    print(f"  Empty Weight          We   : {res['W_empty']:>10.1f} lbf")
    print(f"  Fuel Weight           Wf   : {res['W_fuel']:>10.1f} lbf")
    print(f"  Payload               Wp   : {res['W_payload']:>10.1f} lbf")
    print(f"  We/W0                      : {res['W_empty']/res['W0']:>10.4f}")
    print(f"  Wf/W0                      : {res['W_fuel']/res['W0']:>10.4f}")

    print(f"\n  Aerodynamics")
    print(f"  (L/D)_max                  : {res['LD_max']:>10.2f}")
    print(f"  (L/D) at loiter            : {res['loiter_LD']:>10.2f}")

    print(f"\n  Mission")
    print(f"  Climb covers               : {res['x_climb_nmi']:>10.1f} nmi")
    print(f"  Cruise range (net)         : {res['R_cruise_actual']:>10.1f} nmi")

    print(f"\n  Fuel Breakdown")
    print(f"  {'-'*44}")
    for phase, wt in res['fuel_breakdown'].items():
        pct = 100 * wt / res['W_fuel']
        print(f"  {phase:<20} : {wt:>8.1f} lbf  ({pct:5.1f}%)")
    print(f"  {'Total':<20} : {res['W_fuel']:>8.1f} lbf")

    print(f"\n  Component Weight Breakdown")
    print(f"  {'-'*54}")
    print(f"  {'Group / Component':<28}  {'Weight':>8}   {'% We':>6}")
    print(f"  {'-'*54}")

    def row(label, w):
        pct = 100 * w / res['W_empty']
        print(f"  {label:<28}  {w:>8.1f}   {pct:>5.1f}%")

    print(f"  {'--- Structures ---':<28}  {cw['W_structures']:>8.1f}")
    row("  Wing",          cw['W_wing'])
    row("  Vertical Tail", cw['W_vt'])
    row("  Fuselage",      cw['W_fuselage'])
    row("  AIS",           cw['W_ais'])
    row("  Main LG",       cw['W_mlg'])
    row("  Nose LG",       cw['W_nlg'])
    row("  Eng Section",   cw['W_ensect'])
    row("  Canopy",        cw['W_canopy'])
    row("  Arr Gear",      cw['W_arrgear'])
    row("  Cat Gear",      cw['W_catgear'])
    row("  Handling Gear", cw['W_handgear'])

    print(f"  {'--- Propulsion ---':<28}  {cw['W_propulsion']:>8.1f}")
    row("  Engine (uninstalled)", cw['W_en'])
    row("  Engine Mounts",  cw['W_emounts'])
    row("  Engine Controls",cw['W_encontrols'])
    row("  Engine Cooling", cw['W_encool'])
    row("  Oil Cooling",    cw['W_oilcool'])
    row("  Tailpipe",       cw['W_tailpipe'])
    row("  Starter",        cw['W_starter'])
    row("  Fuel Sys/Tanks", cw['W_fuelsystanks'])

    print(f"  {'--- Equipment ---':<28}  {cw['W_equipment']:>8.1f}")
    row("  Avionics",       cw['W_avionics'])
    row("  Flight Controls",cw['W_flightcont'])
    row("  Instruments",    cw['W_instr'])
    row("  Hydraulics",     cw['W_hydraulics'])
    row("  Electrical",     cw['W_elec'])
    row("  AC/Anti-ice",    cw['W_acai'])
    row("  Furnishings",    cw['W_furn'])
    row("  Firewall",       cw['W_firewall'])

    print(f"  {'-'*54}")
    print(f"  {'TOTAL EMPTY':<28}  {res['W_empty']:>8.1f}")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# 14.  PLOTTING
# ---------------------------------------------------------------------------

def plot_results(res: dict, p: AircraftParams, label: str = ''):
    """
    Plot weight fraction and cumulative fuel burn vs range over the full
    flight regime (climb + cruise), plus a fuel breakdown pie chart.
    """
    ch = res['climb_history']
    cr = res['cruise_history']
    W0 = res['W0']

    # ------------------------------------------------------------------
    # Build unified range axis: climb phase then cruise phase
    # ------------------------------------------------------------------
    # Climb: x_nmi is already cumulative from brake-release
    climb_R   = ch['x_nmi']                          # nmi
    climb_W   = ch['W']                              # lbf at end of each step
    climb_fuel = ch['fuel_cum']                       # cumulative fuel burned [lbf]

    # Cruise: R_nmi is cumulative within the cruise segment; offset by climb distance
    cruise_offset = climb_R[-1] if len(climb_R) > 0 else 0.0
    cruise_R  = cr['R_nmi'] + cruise_offset           # nmi (offset to follow climb)
    cruise_W  = cr['W']                              # lbf
    cruise_fuel = cr['fuel_cum'] + climb_fuel[-1]    # continue cumulative from climb end

    # Concatenate for full-mission arrays
    full_R    = np.concatenate([climb_R,  cruise_R])
    full_W    = np.concatenate([climb_W,  cruise_W])
    full_fuel = np.concatenate([climb_fuel, cruise_fuel])

    # Weight fraction W/W0
    full_Wfrac = full_W / W0

    # Phase boundary for shading
    climb_end_R = cruise_offset

    # ------------------------------------------------------------------
    # Figure: 1×3  (weight fraction | fuel burn | pie)
    # ------------------------------------------------------------------
    title_suffix = f' — {label}' if label else ''
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(f'Refined Sizing — Mission Profile{title_suffix}',
                 fontsize=13, fontweight='bold')

    # ---- Plot 1: Weight fraction vs range ----
    ax = axes[0]
    ax.axvspan(0, climb_end_R, alpha=0.08, color='steelblue', label='Climb')
    ax.axvspan(climb_end_R, full_R[-1], alpha=0.08, color='darkorange', label='Cruise')
    ax.plot(full_R, full_Wfrac, 'b-', linewidth=2)
    ax.set_xlabel('Cumulative Range [nmi]')
    ax.set_ylabel('Weight Fraction  W / W₀')
    ax.set_title('Weight Fraction vs Range')
    ax.set_ylim(0.5, 1.05)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ---- Plot 2: Cumulative fuel burned vs range ----
    ax = axes[1]
    ax.axvspan(0, climb_end_R, alpha=0.08, color='steelblue', label='Climb')
    ax.axvspan(climb_end_R, full_R[-1], alpha=0.08, color='darkorange', label='Cruise')
    ax.plot(full_R, full_fuel / 1000, 'r-', linewidth=2)
    ax.set_xlabel('Cumulative Range [nmi]')
    ax.set_ylabel('Fuel Burned [klbf]')
    ax.set_title('Cumulative Fuel Consumption vs Range')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ---- Plot 3: Fuel breakdown pie chart ----
    ax = axes[2]
    breakdown = res['fuel_breakdown']
    labels = list(breakdown.keys())
    values = list(breakdown.values())
    labels_f = [l for l, v in zip(labels, values) if v > 1.0]
    values_f = [v for v in values if v > 1.0]
    wedges, _ = ax.pie(values_f, labels=None, startangle=90)
    legend_labels = [f'{l}  {100*v/sum(values_f):.1f}%' for l, v in zip(labels_f, values_f)]
    ax.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(0.95, 0.5),
              fontsize=8, frameon=True)
    ax.set_title('Fuel Breakdown by Phase')

    plt.tight_layout()
    plt.show(block=False)


# ---------------------------------------------------------------------------
# 15.  TRADE STUDY HELPER
# ---------------------------------------------------------------------------

def trade_study_AR(p: AircraftParams,
                   AR_range=None,
                   W0_init: float = 45000.0) -> None:
    """
    Run the sizing loop for a range of aspect ratios and plot W0 vs AR.
    This demonstrates the slide 11-12 point about needing physics-based
    wing weight to see the real tradeoff between aero and structures.
    """
    if AR_range is None:
        AR_range = np.linspace(3.0, 10.0, 8)

    W0_results = []
    DERIVED = {'k', 'TW_target'}   # fields computed in __post_init__, not init args
    for AR_val in AR_range:
        # Pass AR_val directly into the constructor so __post_init__ recomputes k
        p_copy = AircraftParams(**{**{f.name: getattr(p, f.name)
                                      for f in p.__dataclass_fields__.values()
                                      if f.name not in DERIVED},
                                   'AR': AR_val})
        res = run_sizing(p_copy, W0_guess=W0_init, verbose=False)
        W0_results.append(res['W0'])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(AR_range, np.array(W0_results) / 1000, 'bo-', linewidth=2, markersize=6)
    ax.set_xlabel('Aspect Ratio (AR)', fontsize=12)
    ax.set_ylabel('Gross Takeoff Weight W0 [klbf]', fontsize=12)
    ax.set_title('Trade Study: AR vs W0\n(Physics-based wing weight)', fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show(block=False)


# ---------------------------------------------------------------------------
# 16.  ENTRY POINT
# ---------------------------------------------------------------------------

if __name__ == '__main__':

    # --- Payload component weights ---
    W_AIM120C  = 358    # lb (from: https://www.navair.navy.mil/product/AMRAAM)
    W_AIM9X    = 186    # lb (from: https://www.navair.navy.mil/product/AIM-9X-Sidewinder)
    W_JDAM     = 1015   # lb (from: https://en.wikipedia.org/wiki/Mark_83_bomb)
    W_avionics = 2500   # lb (avionics/sensors)

    num_AIM120C = 6
    num_AIM9X   = 2
    num_JDAM    = 4

    W_payload_a2a    = num_AIM120C * W_AIM120C + num_AIM9X * W_AIM9X + W_avionics  # lb
    W_payload_strike = num_JDAM * W_JDAM + num_AIM9X * W_AIM9X + W_avionics        # lb

    print(f"W_payload_air2air: {W_payload_a2a} lb")
    print(f"W_payload_strike:  {W_payload_strike} lb")

    # --- Set up your aircraft parameters here ---
    # These are placeholder values representative of an F/A-18C.
    # Replace with your team's PDR values.
    params = AircraftParams(
        # --- Geometry / aero ---
        S_ref      = 600.0,
        AR         = 2.09,
        e          = 0.82,
        tc         = 0.05,
        taper      = 0.30,
        sweep      = 40.0,
        CD0        = 0.0165,
        # --- Propulsion ---
        T_sl_max   = 43000.0,
        n_engines  = 1,
        ct_cruise  = 0.866,
        ct_loiter  = 0.80,
        ct_sl      = 0.40,
        idle_frac  = 0.05,
        # --- Payloads ---
        W_payload_a2a    = W_payload_a2a,
        W_payload_strike = W_payload_strike,
        # --- Mission ---
        h_cruise     = 40000.0,
        M_cruise     = 0.85,
        R_cruise     = 700 * 2,
        n_cruise_seg = 20,
        E_loiter     = 0.3333,
        # --- Component weight parameters (from Component_Weights_and_CG.py) ---
        K_cb   = 1.0,
        K_d    = 2.75,
        K_dw   = 0.768,
        K_dwf  = 0.774,
        K_mc   = 1.45,
        K_rht  = 1.047,
        K_tpg  = 1.0,
        K_vg   = 1.62,
        K_vs   = 1.0,
        K_vsh  = 1.0,
        L_d    = 5.809,
        L_ec   = 10.503,
        L_sh   = 220/12,
        L_t    = 8.72490,
        L_m    = 60.0,
        L_n    = 60 + (2.704*12),
        S_fw   = 400.0,
        S_vt   = 100.0,
        S_r    = 100.0,
        S_csw  = 200.0,
        H_tH_v = 0.0,
        L_fus  = 40.28139,
        D_fus  = 5.4,
        W_fus  = 13.232,
        D_e    = 43/12,
        V_pr   = 120.0,
        AR_vt    = 1.00176,
        tc_vt    = 0.100,
        taper_vt = 0.25,
        sweep_vt = 45.0,
        V_i    = 120 * 7.48052,
        N_t    = 1.0,
        SFC    = 0.886,
        W_en   = 6422.0,
        N_c    = 1.0,
        N_ci   = 1.0,
        N_s    = 1.0,
        N_u    = 10.0,
        N_nw   = 2.0,
        N_gear = 3.0,
        R_kva  = 160.0,
        W_uav  = 2500.0,
        wing_fudge     = 0.85,
        tail_fudge     = 0.83,
        fuselage_fudge = 0.90,
        AIS_fudge      = 0.85,
        Nz       = 10.5,
        limloadf = 7.0,
        W_dg_PDR = 43200.0,   # PDR gross weight — anchors T/W ratio
    )

    # --- Run the sizing loop for both loadouts ---
    print("\n--- Air-to-Air Loadout ---")
    results_a2a = run_sizing(params, W_payload=params.W_payload_a2a,
                             W0_guess=45000.0, tol=1.0, verbose=True)
    print_results(results_a2a, params)

    print("\n--- Ground Strike Loadout ---")
    results_strike = run_sizing(params, W_payload=params.W_payload_strike,
                                W0_guess=45000.0, tol=1.0, verbose=True)
    print_results(results_strike, params)

    # --- Generate plots for both loadouts ---
    plot_results(results_a2a,    params, label='Air-to-Air')
    plot_results(results_strike, params, label='Ground Strike')

    # --- AR trade study ---
    print("\nRunning AR trade study...")
    trade_study_AR(params, AR_range=np.linspace(3.0, 9.0, 10), W0_init=45000.0)

    # Keep all figure windows open until manually closed
    plt.show()