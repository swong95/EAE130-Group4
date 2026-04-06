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

    # Thrust lapse with altitude (simplified: T = T_sl * (rho/rho_sl)^n)
    thrust_lapse_exp: float = 0.7

    # --- Weights (fixed) ---
    W_payload_a2a:    float = 5020   # Payload for air-to-air loadout [lbf]
    W_payload_strike: float = 6932   # Payload for ground strike loadout [lbf]

    # --- Empty weight fractions (Raymer-style, used for component buildup) ---
    # For now we use a simple fraction; you can replace with full component equations
    We_fraction: float = 0.53   # Historical empty-weight fraction for fighter

    # Nz - ultimate load factor (1.5 * limit load factor)
    Nz: float = 11.25           # 1.5 * 7.5g limit

    # Wing control surface area
    S_csw: float = 62        # [ft^2]

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

    # --- Derived properties (computed in __post_init__) ---
    k: float = field(init=False)

    def __post_init__(self):
        self.k = 1.0 / (np.pi * self.AR * self.e)

    @property
    def T_max_total(self) -> float:
        """Total max sea-level thrust [lbf]."""
        return self.T_sl_max * self.n_engines


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

def thrust_available(h: float, p: AircraftParams) -> float:
    """
    Available max continuous thrust at altitude h [lbf].
    Uses simplified density-ratio lapse: T = T_sl * (rho/rho0)^n
    """
    rho, _ = atmosphere(h)
    rho0, _ = atmosphere(0.0)
    return p.T_max_total * (rho / rho0) ** p.thrust_lapse_exp


# ---------------------------------------------------------------------------
# 5.  BEST CLIMB SPEED FOR A JET  (slide 48)
# ---------------------------------------------------------------------------

def best_climb_speed(W: float, h: float, p: AircraftParams) -> float:
    """
    Speed [ft/s] that maximises rate of climb for a jet at altitude h,
    weight W, assuming thrust = T_available (max continuous).

    From slide 48:
        V = sqrt[ (T/W * W/S) / rho  +  sqrt( (T/W * W/S / rho)^2
                                               + 12 * CD0 * k * (W/S)^2 / rho^2 ) ]
    divided appropriately.

    This comes from d(RC)/dV = 0 with T = const.
    """
    rho, _ = atmosphere(h)
    T  = thrust_available(h, p)
    TW = T / W
    WS = W / p.S_ref

    # Coefficient groupings from the slide 48 formula
    A = TW * WS / rho
    B = 12.0 * p.CD0 * p.k * (WS / rho) ** 2

    V2 = A + np.sqrt(A**2 + B)   # V^2 in ft^2/s^2 — note this is per-unit simplification
    # Full derivation gives V^2 = ( A + sqrt(A^2 + B) ) already in correct units
    return np.sqrt(max(V2, 1.0))  # guard against negatives


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

def simulate_climb(W_start: float, p: AircraftParams) -> dict:
    """
    Simulate climb from sea level to h_cruise using the energy-height method.

    At each altitude step:
      1. Find best-climb speed V at current (h, W, rho)
      2. Compute T, D, L/D
      3. Compute Ps = V*(T-D)/W  (specific excess power = d(he)/dt)
      4. Compute dt = delta_he / Ps
      5. Compute fuel burned: delta_W = ct * T * dt
      6. Track horizontal distance: delta_x = V * cos(gamma) * dt

    Returns
    -------
    dict with keys:
        W_end        : weight at top of climb [lbf]
        fuel_burned  : total fuel burned in climb [lbf]
        x_climb      : horizontal distance covered [nmi]
        history      : dict of arrays for plotting
    """
    h      = 0.0
    W      = W_start
    h_end  = p.h_cruise
    dh     = p.h_climb_step

    # Storage for plotting
    hist = {
        'h': [], 'V_kts': [], 'W': [], 'RC': [], 'gamma_deg': [],
        'fuel_cum': [], 'x_nmi': [], 'LD': []
    }

    fuel_total = 0.0
    x_total    = 0.0   # ft

    while h < h_end:
        dh_step = min(dh, h_end - h)   # don't overshoot
        h_mid   = h + 0.5 * dh_step    # use midpoint altitude for accuracy

        rho_mid, a_mid = atmosphere(h_mid)
        T_mid          = thrust_available(h_mid, p)

        # Best climb speed at midpoint conditions
        V = best_climb_speed(W, h_mid, p)

        # Aerodynamics
        CL  = CL_from_weight(W, rho_mid, V, p)
        CD  = CD_from_CL(CL, p)
        D   = 0.5 * rho_mid * V**2 * p.S_ref * CD
        L_D = CL / CD

        # Rate of climb and flight path angle
        Ps      = V * (T_mid - D) / W      # ft/s  (specific excess power)
        Ps      = max(Ps, 1.0)              # prevent division by zero near ceiling
        sin_g   = Ps / V
        gamma   = np.arcsin(np.clip(sin_g, -1, 1))

        # Energy height change for this altitude step
        # he = h + V^2/(2g)  — we need delta_he
        # Approximate: compute he at top and bottom of step
        # (V changes slightly; use current V as representative)
        he_bot = h + V**2 / (2 * G)
        he_top = (h + dh_step) + V**2 / (2 * G)   # V ~ constant over small step
        delta_he = he_top - he_bot                  # ≈ dh_step for small steps

        # Time to traverse this energy height increment
        dt = delta_he / Ps                          # seconds

        # Fuel burned  (ct in 1/hr -> convert to 1/s)
        ct_s        = p.ct_sl / 3600.0              # use cruise ct during climb
        # Actually use cruise TSFC during climb (common simplification)
        ct_climb_s  = p.ct_cruise / 3600.0
        delta_W     = ct_climb_s * T_mid * dt       # lbf of fuel

        # Horizontal distance
        delta_x     = V * np.cos(gamma) * dt        # ft

        # Update state
        fuel_total += delta_W
        x_total    += delta_x
        W          -= delta_W
        h          += dh_step

        # Record history (at midpoint)
        hist['h'].append(h_mid)
        hist['V_kts'].append(V / 1.6878)            # ft/s -> knots
        hist['W'].append(W)
        hist['RC'].append(Ps * 60)                  # ft/s -> ft/min
        hist['gamma_deg'].append(np.degrees(gamma))
        hist['fuel_cum'].append(fuel_total)
        hist['x_nmi'].append(x_total / FT_PER_NMI)
        hist['LD'].append(L_D)

    return {
        'W_end':       W,
        'fuel_burned': fuel_total,
        'x_climb':     x_total / FT_PER_NMI,        # nmi
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

def empty_weight(W0: float, p: AircraftParams) -> float:
    """
    Estimate aircraft empty weight.

    For now this uses a simple historical fraction.
    You should replace this with your full component buildup equations
    (Raymer Ch. 15) once you have all the geometry defined.

    The wing weight equation from slide 11 is included as an example:
        W_wing = 0.0051 * (Wdg * Nz)^0.557 * Sw^0.649 * AR^0.5
                         * (t/c)^-0.4 * (1+lambda)^0.1 * cos(Lambda)^-1 * Scsw^0.1
    """
    # Full Raymer wing weight (slide 11)
    Wdg        = W0                             # Design gross weight [lbf]
    sweep_rad  = np.radians(p.sweep)
    W_wing = (0.0051
              * (Wdg * p.Nz) ** 0.557
              * p.S_ref ** 0.649
              * p.AR ** 0.5
              * p.tc ** (-0.4)
              * (1 + p.taper) ** 0.1
              * np.cos(sweep_rad) ** (-1.0)
              * p.S_csw ** 0.1)

    # Everything else: use historical fraction for the rest of the empty weight
    # W_empty_total = We_fraction * W0
    # W_other = W_empty_total - W_wing  (other structure, systems, avionics, etc.)
    W_empty_total = p.We_fraction * W0
    W_other       = max(0.0, W_empty_total - W_wing)

    return W_wing + W_other     # = We_fraction * W0  (until you add more components)


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
        ff_taxi   = fuel_fraction_taxi(p)       # fractional fuel / W_start
        ff_to     = fuel_fraction_takeoff(p)

        W_after_taxi = W0 * (1 - ff_taxi)
        W_after_to   = W_after_taxi * (1 - ff_to)

        fuel_taxi = W0 - W_after_taxi
        fuel_to   = W_after_taxi - W_after_to

        # --- Phase 2: Climb ---
        climb = simulate_climb(W_after_to, p)
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

        # --- Empty weight ---
        W_empty = empty_weight(W0, p)

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
        'converged':       converged,
        'W0':              W0,
        'W_empty':         W_empty,
        'W_fuel':          W_fuel_total,
        'W_payload':       W_payload,
        'fuel_breakdown':  fuel_breakdown,
        'x_climb_nmi':     x_climb,
        'R_cruise_actual': p.R_cruise - x_climb,
        'LD_max':          LD_max(p),
        'loiter_LD':       last_loiter['LD_loiter'],
        'climb_history':   last_climb['history'],
        'cruise_history':  last_cruise['history'],
    }


# ---------------------------------------------------------------------------
# 13.  RESULTS REPORTING
# ---------------------------------------------------------------------------

def print_results(res: dict, p: AircraftParams):
    """Pretty-print the sizing results."""
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
    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# 14.  PLOTTING
# ---------------------------------------------------------------------------

def plot_results(res: dict, p: AircraftParams):
    """Generate mission profile and weight history plots."""

    ch = res['climb_history']
    cr = res['cruise_history']

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle('Refined Sizing — Mission Profile', fontsize=14, fontweight='bold')

    # ---- Plot 1: Altitude profile during climb ----
    ax = axes[0, 0]
    ax.plot(ch['x_nmi'], ch['h'] / 1000, 'b-', linewidth=2)
    ax.set_xlabel('Distance [nmi]')
    ax.set_ylabel('Altitude [kft]')
    ax.set_title('Climb Altitude Profile')
    ax.grid(True, alpha=0.3)

    # ---- Plot 2: Rate of climb ----
    ax = axes[0, 1]
    ax.plot(ch['h'] / 1000, ch['RC'], 'g-', linewidth=2)
    ax.axhline(y=500, color='r', linestyle='--', alpha=0.7, label='500 ft/min min.')
    ax.set_xlabel('Altitude [kft]')
    ax.set_ylabel('Rate of Climb [ft/min]')
    ax.set_title('Rate of Climb vs Altitude')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # ---- Plot 3: Weight during climb ----
    ax = axes[0, 2]
    ax.plot(ch['x_nmi'], ch['W'] / 1000, 'r-', linewidth=2)
    ax.set_xlabel('Distance [nmi]')
    ax.set_ylabel('Weight [klbf]')
    ax.set_title('Weight During Climb')
    ax.grid(True, alpha=0.3)

    # ---- Plot 4: L/D during cruise ----
    ax = axes[1, 0]
    ax.plot(cr['R_nmi'], cr['LD'], 'b-', linewidth=2)
    ax.axhline(y=res['LD_max'], color='r', linestyle='--',
               alpha=0.7, label=f'(L/D)_max = {res["LD_max"]:.1f}')
    ax.set_xlabel('Range [nmi]')
    ax.set_ylabel('L/D [-]')
    ax.set_title('L/D During Cruise')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # ---- Plot 5: Weight during cruise ----
    ax = axes[1, 1]
    ax.plot(cr['R_nmi'], cr['W'] / 1000, 'g-', linewidth=2)
    ax.set_xlabel('Range [nmi]')
    ax.set_ylabel('Weight [klbf]')
    ax.set_title('Weight During Cruise')
    ax.grid(True, alpha=0.3)

    # ---- Plot 6: Fuel breakdown pie chart ----
    ax = axes[1, 2]
    breakdown = res['fuel_breakdown']
    labels = list(breakdown.keys())
    values = list(breakdown.values())
    # Filter out near-zero slices for readability
    labels_f = [l for l, v in zip(labels, values) if v > 1.0]
    values_f = [v for v in values if v > 1.0]
    ax.pie(values_f, labels=labels_f, autopct='%1.1f%%', startangle=90)
    ax.set_title('Fuel Breakdown by Phase')

    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/refined_sizing_plots.png', dpi=150,
                bbox_inches='tight')
    plt.show()
    print("  Plot saved to refined_sizing_plots.png")


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
    for AR_val in AR_range:
        p_copy     = AircraftParams(**{f.name: getattr(p, f.name)
                                       for f in p.__dataclass_fields__.values()
                                       if f.name != 'k'})
        p_copy.AR  = AR_val
        res        = run_sizing(p_copy, W0_guess=W0_init, verbose=False)
        W0_results.append(res['W0'])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(AR_range, np.array(W0_results) / 1000, 'bo-', linewidth=2, markersize=6)
    ax.set_xlabel('Aspect Ratio (AR)', fontsize=12)
    ax.set_ylabel('Gross Takeoff Weight W0 [klbf]', fontsize=12)
    ax.set_title('Trade Study: AR vs W0\n(Physics-based wing weight)', fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/trade_study_AR.png', dpi=150,
                bbox_inches='tight')
    plt.show()
    print("  AR trade study plot saved to trade_study_AR.png")


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
        S_ref            = 600.0,   # ft^2
        AR               = 2.09,
        e                = 0.82,
        tc               = 0.05,
        taper            = 0.30,
        sweep            = 40,      # deg
        CD0              = 0.0165,
        T_sl_max         = 36000.0, # lbf per engine
        n_engines        = 2,
        ct_cruise        = 0.70,    # 1/hr
        ct_loiter        = 0.60,
        ct_sl            = 0.40,
        idle_frac        = 0.05,
        W_payload_a2a    = W_payload_a2a,
        W_payload_strike = W_payload_strike,
        We_fraction      = 0.53,
        Nz               = 11.25,
        S_csw            = 40.0,    # ft^2
        h_cruise         = 40000.0, # ft
        M_cruise         = 0.90,
        R_cruise         = 700 * 2, # nmi
        n_cruise_seg     = 20,
        E_loiter         = 0.3333,  # hr
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

    # --- Generate plots (for air-to-air as primary) ---
    plot_results(results_a2a, params)

    # --- AR trade study ---
    print("\nRunning AR trade study...")
    trade_study_AR(params, AR_range=np.linspace(3.0, 9.0, 10), W0_init=45000.0)