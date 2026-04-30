import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# weight data
MTOW = 48891.4578 # 52027.7578 lbf for strike / 48891.4578 lbf for a2a from jonas' component_weights_and_cg
OEW  = 26614.557835667372   # lbf operating empty weight (empty weight from ^)
MFW = 18967.3 # max fuel weight capacity W_fuel_a2a = 18967.3 lbf W_fuel_strike = 20191.6 
PW = 2520 #Payload A2A Weight: 2520 Payload A2A Strike: 4432
MZFW = MTOW - MFW  # lb max zero fuel weight (check if this is correct?)
MPW = MZFW - OEW 
RF = 1800   #lbf: find a good source for this? https://info.publicintelligence.net/F18-EF-200.pdf

# flight params from A2 improved drag polars
rho_c  = 0.0008893             # slug/ft^3 (air density at 30,000 ft.)
S      = 600     # ft2       wing reference area
CL     = 0.36703409204468285       #           cruise lift coefficient
LD     = 0.36703409204468285/0.04378150000000001     #           cruise L/D  (lift-to-drag ratio)
CD     = CL / LD    #           cruise drag coefficient  (= 0.026)
TSFC   = 0.8 # specific fuel consumption (Placeholder? from some random website for FW-135)
ct     = TSFC / 3600            # lb/lb/s   (convert for the Breguet integral)

# breuget range function
def breguet(W0, W1):
    """
    Jet Breguet range (constant altitude, constant CL).
    R = 2*sqrt(2/rho*S) * (1/ct) * (CL^0.5/CD) * (sqrt(W0) - sqrt(W1))
    W0, W1 in lb  ->  returns range in nautical miles
    """
    # FT to NM conversion factor
    FT_PER_NM = 6076.12  # ft / nautical mile

    # Compute range in feet using the Breguet formula
    R_ft = (2 * np.sqrt(2 / (rho_c * S))
              * (1 / ct)
              * (CL**0.5 / CD)
              * (W0**0.5 - W1**0.5))
    range_nm = R_ft / FT_PER_NM # convert ft to NM

    return range_nm

payload_A = MPW # point A has the max payload weight
range_A   = 0

print(f'Point A')
print(f'  Payload = {payload_A:,} lb  (= MPW)')
print(f'  Range   = {range_A} nm')

FW_B      = MTOW - OEW - MPW      # fuel that fits after max payload, at MTOW
W_res_B   = RF          # reserve (not burned for range)
payload_B = MPW
W0_B      = MTOW
W1_B      = OEW + MPW + W_res_B
range_B   = breguet(W0_B, W1_B)

print(f'Point B')
print(f'  Fuel loaded  = MTOW - OEW - MPW = {MTOW} - {OEW} - {MPW} = {FW_B:,} lb')
# print(f'  Reserve      = {RFF:.0%} x {FW_B:,} = {W_res_B:,.0f} lb')
print(f'  W0 = MTOW                       = {W0_B:,} lb')
print(f'  W1 = OEW + MPW + Wres           = {OEW} + {MPW} + {W_res_B:.0f} = {W1_B:,.0f} lb')
print(f'  Payload = {payload_B:,} lb')
print(f'  Range   = {range_B:,.0f} nm')

FW_C      = MFW
payload_C = MTOW - OEW - MFW      # payload reduced to stay within MTOW
W_res_C   = RF
W0_C      = MTOW
W1_C      = OEW + payload_C + W_res_C
range_C   = breguet(W0_C, W1_C)

print(f'Point C')
print(f'  Fuel loaded  = MFW             = {FW_C:,.0f} lb')
print(f'  Payload_C    = MTOW-OEW-MFW    = {MTOW}-{OEW}-{FW_C:.0f} = {payload_C:,.0f} lb  (reduced!)')
# print(f'  Reserve      = {RFF:.0%} x {FW_C:,.0f} = {W_res_C:,.0f} lb')
print(f'  W0 = MTOW                      = {W0_C:,} lb')
print(f'  W1 = OEW + Payload_C + Wres    = {OEW} + {payload_C:.0f} + {W_res_C:.0f} = {W1_C:,.0f} lb')
print(f'  Payload = {payload_C:,.0f} lb')
print(f'  Range   = {range_C:,.0f} nm')

FW_D      = MFW
payload_D = 0
W_res_D   = RF
W0_D      = OEW + FW_D            # below MTOW — no payload means lighter aircraft
W1_D      = OEW + W_res_D
range_D   = breguet(W0_D, W1_D)

print(f'Point D')
print(f'  Fuel loaded  = MFW             = {FW_D:,.0f} lb')
print(f'  W0 = OEW + MFW                 = {OEW} + {FW_D:.0f} = {W0_D:,.0f} lb  (< MTOW = {MTOW:,})')
# print(f'  Reserve      = {RFF:.0%} x {FW_D:,.0f} = {W_res_D:,.0f} lb')
print(f'  W1 = OEW + Wres                = {OEW} + {W_res_D:.0f} = {W1_D:,.0f} lb')
print(f'  Payload = {payload_D} lb')
print(f'  Range   = {range_D:,.0f} nm')


target_range_nm = 1400
FT_PER_NM = 6076.12

# Rearrange Breguet to find required W1 for a given Range
term_a = 2 * np.sqrt(2 / (rho_c * S)) * (1 / ct) * (CL**0.5 / CD)
W0_T = MTOW
W1_T_required = (np.sqrt(W0_T) - (target_range_nm * FT_PER_NM) / term_a)**2

# Payload = W1 - OEW - Reserve
payload_T = PW
range_T = target_range_nm 
FW_T = MTOW - OEW - payload_T # Fuel loaded for mission T

print(f'Point T (Fixed Range)')
print(f'  Target Range = {range_T} nm')
print(f'  Required Payload = {payload_T:,.0f} lb')


summary = pd.DataFrame({
    'Point': ['A', 'B', 'C', 'D', 'T'],
    'Description': [
        'Zero range / max payload',
        'Harmonic range (max payload @ MTOW)',
        'Max-fuel range (full tanks @ MTOW)',
        'Ferry range (zero payload)',
        'Design / target mission',
    ],
    'Payload (lb)': [int(payload_A), int(payload_B), int(round(payload_C)),
                     int(payload_D), int(payload_T)],
    'TOW (lb)':     [int(OEW+payload_A), int(W0_B), int(W0_C),
                     int(round(W0_D)), int(W0_T)],
    'Fuel (lb)':    [0, int(FW_B), int(round(FW_C)), int(round(FW_D)), int(FW_T)],
    'W0 (lb)':      [int(OEW+payload_A), int(W0_B), int(W0_C),
                     int(round(W0_D)), int(W0_T)],
    'W1 (lb)':      [int(OEW+payload_A), int(round(W1_B)), int(round(W1_C)),
                     int(round(W1_D)), int(round(W1_T_required))],
    'Range (nm)':   [int(range_A), int(round(range_B)), int(round(range_C)),
                     int(round(range_D)), int(round(range_T))],
})

fig, ax = plt.subplots(figsize=(12, 7))

# Envelope A -> B -> C -> D
env_r = [range_A,   range_B,   range_C,   range_D]
env_p = [payload_A, payload_B, payload_C, payload_D]
ax.plot(env_r, env_p, color='steelblue', linewidth=2.5,
        marker='s', markersize=8, label='Payload-Range Envelope')

# Design point T
ax.plot(range_T, payload_T, color='crimson', marker='D',
        markersize=10, linestyle='None', label='Design Mission (A2A Payload and 1,400 nm)')
ax.axhline(payload_T, color='crimson', linestyle='--', linewidth=1.1, alpha=0.5)
ax.axvline(range_T,   color='crimson', linestyle='--', linewidth=1.1, alpha=0.5)

# Annotate envelope points
offsets = {'A': (15, 100), 'B': (15, 100),
           'C': (30, -600), 'D': (15, 100)}
for pt, r, p in zip(['A','B','C','D'], env_r, env_p):
    dx, dy = offsets[pt]
    ax.annotate(f'{pt}\n({int(round(r)):,} nm, {int(round(p)):,} lb)',
                xy=(r, p), xytext=(r+dx, p+dy), fontsize=9, va='bottom')

# Annotate Point T and MPW
ax.annotate(f'T: {int(payload_T):,} lb\nat {range_T} nm',
            xy=(range_T, payload_T), xytext=(range_T-200, payload_T+25),
            fontsize=10, color='crimson', fontweight='bold')


ax.set_xlabel('Range (nm)', fontsize=12)
ax.set_ylabel('Payload (lb)', fontsize=12)
ax.set_title('Payload-Range Chart — F/A XX (Air to Air (A2A)) Config)', fontsize=14)
ax.set_xlim(0, max(env_r)+200)
ax.set_ylim(0, MPW+1000)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(fontsize=10, loc='upper right')
plt.tight_layout()
plt.show()

print(summary)