import numpy as np
import matplotlib.pyplot as plt

# Parameters from Raymer Chapter 15 and OPENVSP Model
# (R) = Raymer
# (VSP) = from model
K_cb = 1.0 # main landing gear crossbeam value
K_d = 2.75 # duct constant (R)
K_dw = 0.768 #for delta wing (R)
K_dwf = 0.774 # for delta wing (R)
K_mc = 1.45 # electronics mission completion value (R)
K_rht = 1.047 # for rolling tail (R)(2 horizontal tail surfaces move independently)
K_tpg = 1.0 # tripod gear value
K_vg = 1.62 # for variable geometry of air induction system
K_vs = 1.0 # for variable sweep wing value(R)
K_vsh = 1.0 # for variable sweep wing hydraulics value(R)
L_d = 5.809 # ft length of duct (VSP)
L_ec = 10.503 # ft length from engine front to cockpit (VSP)
L_a = L_ec # ft electrical routing distance, generators to avionics to cockpit
L_m = 60 # in main landing gear length in inches (Placeholder)(not modeled)
L_n = 60 + (2.704*12) # in nose landing gear length in inches (Placeholder)(not modeled)
L_s = 5.809 # ft single duct lenth 
L_sh = 220 / 12 #ft shroud length (Placeholder)(Guess)
L_t = 8.72490 # length tail ft
H_tH_v = 0 # no horizontal tail


W_dg = 51908.0 # Design Weight (From A1_KevTestCode.py)
limloadf = 7 # g
N_z = 1.5 * limloadf # Ultimate load factor (R)
S_w = 600 # ft^2 Trapezoidal Wing Area
S_vt = 100 # vertical tail area ft^2
S_fw = 400 # ft^2 firewall surface area (Placeholder)(not modeled)
S_csw = 200 # Wing control surface area ft^2 (Placeholder need to run in VSPAERO)
S_r = 100 # ft^2 Rudder area
S_cs = S_r + S_csw # ft^2 total area of control surfaces
SFC = 0.886 # specific fuel consumption (Placeholder? from some random website for FW-135)
T = 43000 #lb total engine thrust
T_e = T # lb thrust per engine
AR = 2.08815 # Aspect Ratio
AR_vt = 1.00176 # AR vertical tail
tc_root = 0.125 # root thickness/chord
tc_root_vt = 0.10000 # root t/c for vertical tail
M = 2 # Mach Number 
MAC = 20.84049 #ft
N_c = 1 # number of crew
N_ci = 1.0 # 1.0 single pilot, 1.2 pilot + backseater, 2.0 pilot & copassenger
N_en = 1 # number of engines
N_gen = N_en # number of generators (typically same as num engines)
N_gear = 3 # Assuming this is number of landing gear in total
N_l = N_gear * 1.5 # ultimate landg load factor
N_s = 1 # number of flight control systems
N_t = 1 # number of fuel tanks (may need updating in future)
N_nw = 2 # number of nose wheels
N_u = 10 # number of hydraulic utility functions (typically 5-15) (R)
c_tip = 3.28876 # ft
c_root = 32.88764 # ft
c_tip_vt = 2.29229
c_root_vt = 9.16916
D = 5.4 # ft fuselage depth (height at max point?)
D_e = 43 / 12# ft enigne diameter (using inlet diameter rn max is 46 in)
L = 40.28139 # ft fuselage structural length
lamb = c_tip/c_root
lamb_vt = c_tip_vt/c_root_vt
biglamb = 40 * np.pi/180 # sweep angle at 25% MAC in radians
biglamb_vt = 45 * np.pi/180 # vertical tail sweep angle at 25% MAC in radians
R_kva = 160 # system electrical rating kv*A (110-160 for fighters & bombers) (R)
V_i = 120* 7.48052 # gal integral tanks volume
V_t = V_i # gal total fuel volume (for later, V_i + other tanks and lines)
V_p = V_i # gal self-sealing "protected" tanks volume
V_pr = 120 # ft^3 volume of pressurized section
W = 13.232 # ft fuselage structural width
W_l = W_dg * 0.78 # design landing gross weight (recheck)
W_en = 6422 # lb engine weight (from Simple Flying and Wikipedia)
W_uav = 2500 # lb uninstalled avionics weight (RFP)

wing_fudge = 0.85 # fudge factor for composites
tail_fudge = 0.83 # fudge factor for composites
fuselage_fudge = 0.90 # fudge factor for composites
AIS_fudge = 0.85 # Air induction system fudge factor for composites

# Wing Weight
W_wing = 0.0103*K_dw*K_vs*((W_dg*N_z)**0.5)*(S_w**0.622)*(AR**0.785)*(tc_root**-0.4)\
    * ((1+lamb)**0.05) * (np.cos(biglamb)**-1.0) * (S_csw**0.04)
W_wing = W_wing*wing_fudge
print("Wing Weight:", W_wing)

# Vertical Tail Weight
W_vt = 0.452*K_rht*((1 + H_tH_v)**0.5) * ((W_dg*N_z)**0.488) * (S_vt**0.718) * (M**0.341)\
    * (L_t**-1.0) * ((1+S_r/S_vt)**0.348) * (AR_vt**0.223) * ((1+lamb_vt)**0.25)\
    * ((np.cos(biglamb_vt))**-0.323)
W_vt = W_vt*tail_fudge
print("Vertical Tail Weight:", W_vt)

# Fuselage Weight
W_fuselage = 0.499*K_dwf*(W_dg**0.35)*(N_z**0.25)*(L**0.5)*(D**0.849)*(W*0.685)
W_fuselage = W_fuselage*fuselage_fudge
print("Fuselage Weight:", W_fuselage)

# Main Landing Gear Weight
W_mlg = K_cb*K_tpg*((W_l*N_l)**0.25)*(L_m**0.973)
print("Main Landing Gear Weight:", W_mlg)

# Nose Landing Gear Weight
W_nlg = ((W_l*N_l)**0.290)*(L_n**0.5)*(N_nw**0.525)
print("Nose Landing Gear Weight:", W_nlg)

# Engine Mount Weights
W_emounts = 0.013*(N_en**0.795)*(T**0.579)*N_z
print("Engine Mount Weight:", W_emounts)

# Firewall Weight
W_firewall = 1.13*S_fw
print("Firewall Weight:", W_firewall)

# Engine Section Weight
W_ensect = 0.01*(W_en**0.717)*N_en*N_z
print("Engine Section Weight:", W_ensect)

# Air Induction System Weight
W_ais = 13.29*K_vg*(L_d**0.643)*(K_d**0.182)*(N_en**1.498)*((L_s/L_d)**-0.373)*D_e
W_ais = W_ais*AIS_fudge
print("Air Induction System Weight:", W_ais)

#Tailpipe Weight
W_tailpipe = 3.5 * D_e*L_sh*N_en
print("Tailpipe Weight:", W_tailpipe)

# Engine Cooling Weight
W_encool = 4.55*D_e*L_sh*N_en
print("Engine Cooling Weight:", W_encool)

# Oil Cooling Weight
W_oilcool = 37.82*(N_en**1.023)
print("Oil Cooling Weight:", W_oilcool)

# Engine Controls Weight
W_encontrols = 10.5*(N_en**1.008)*(L_ec**0.222)
print("Engine Controls Weight:", W_encontrols)

# Starter (Pneumatic) Weight
W_starter = 0.025*(T_e**0.760)*(N_en**0.72)
print("Starter (Pneumatic) Weight:", W_starter)

# Fuel System and Tanks Weight
W_fuelsystanks = 7.45*(V_t**0.47)*((1+(V_i/V_t))**-0.095) * (1+(V_p/V_t))*(N_t**0.066)\
    *(N_en**0.052)*((T*SFC/1000)**0.249)
print("Fuel System and Tanks Weight:", W_fuelsystanks)

# Flight Controls Weight
W_flightcont = 36.28*(M**0.003)*(S_cs**0.489)*(N_s**0.484)*(N_c**0.127)
print("Flight Controls Weight:", W_flightcont)

# Instruments Weight
W_instr = 8.0 + 36.37*(N_en**0.676)*(N_t**0.237) + 26.4*((1+N_ci)**1.356)
print("Instruments Weight:", W_instr)

# Hydraulics Weight
W_hydraulics = 37.23*K_vsh*(N_u**0.664)
print("Hydraulics Weight:", W_hydraulics)

# Electrical Weight
W_elec = 172.2*K_mc*(R_kva**0.152)*(N_c**0.10)*(L_a**0.10)*(N_gen**0.091)
print("Electrical Weight:", W_elec)

# Avionics Weight
W_avionics = 2.117*(W_uav**0.933)
print("Avionics Weight:", W_avionics)

# Furnishings Weight
W_furn = 217.6 * N_c
print("Furnishings Weight:", W_furn)

# Air Conditioning and Anti-ice Weight
W_acai = 201.6*(((W_uav+(200*N_c))/1000)**0.735)
print("Air Conditioning and Anti-ice Weight:", W_acai)

# Handling Gear Weight
W_handgear = 3.2e-4 * W_dg
print("Handling Gear Weight:", W_handgear)

W_canopy = 200 # lb based off F22 rough weight
W_pilot = 200 #lb
W_arrgear = 0.008*W_dg
print("Arresting Gear Weight:", W_arrgear)
W_catgear = 0.003*W_dg
print("Catapult Gear Weight:", W_catgear)

W_AIM120C = 358         # lb (from: https://www.navair.navy.mil/product/AMRAAM)
W_AIM9X  = 186          # lb (from: https://www.navair.navy.mil/product/AIM-9X-Sidewinder)
W_JDAM   = 1015         # lb (from: https://en.wikipedia.org/wiki/Mark_83_bomb)
num_AIM120C = 6
num_AIM9X  = 2
num_JDAM   = 4
W_py_AIM120C = 42.96        
W_py_AIM9X  = 22.32         
W_py_JDAM   = 121.8


# ==================================================================================================================
# CG Calcs
# ==================================================================================================================
origin_offset = 6.4 # ft x direction offset to nose of plane

W_structures = W_wing + W_vt+ W_fuselage+W_ais+W_canopy+W_mlg+W_nlg+W_ensect+W_handgear+W_arrgear+W_catgear
print("Structures Group Weight:", W_structures)

W_propulsion = W_en+W_fuelsystanks+W_emounts+W_encontrols+W_encool+W_oilcool+W_tailpipe+W_starter
print("Propulsion Group Weight:", W_propulsion)

W_payload_a2a = num_AIM120C*W_AIM120C + num_AIM9X*W_AIM9X
W_payload_strike = num_AIM9X*W_AIM9X + num_JDAM*W_JDAM

W_equipment_static = W_avionics+W_firewall+W_flightcont+W_instr+W_hydraulics+W_elec+W_acai+W_furn
print("Equipment Static Group Weight:", W_equipment_static)

W_equipment_a2a = W_payload_a2a + W_equipment_static
print("Equipment A2A Group Weight:", W_equipment_a2a)

W_equipment_strike = W_payload_strike + W_equipment_static
print("Equipment Strike Group Weight:", W_equipment_strike)

# Useful Group Weight Minus Fuel
W_useful = W_pilot + (num_AIM120C*W_py_AIM120C) + (num_AIM9X*W_py_AIM9X) + (num_JDAM*W_py_JDAM)

W_empty = W_structures+W_propulsion+W_equipment_static
print("Empty Weigt:", W_empty)

# W_all_min_fuel = W_structures+W_propulsion+W_equipment_a2a+W_useful
# W_fuel = W_dg-W_all_min_fuel
W_fuel_a2a = 18967.3 #lbf 
W_fuel_strike = 20191.6 #lbf 
# print("Fuel Weight:", W_fuel)

W_useful_report = W_useful + W_fuel_a2a
print("Report Useful Group Weight:", W_useful_report)

W_empty_payload_a2a = W_structures+W_propulsion+W_equipment_static+W_useful+W_fuel_a2a
W_empty_payload_strike = W_structures+W_propulsion+W_equipment_static+W_useful+W_fuel_strike
W_total_a2a = W_structures+W_propulsion+W_equipment_a2a+W_useful+W_fuel_a2a
W_total_strike = W_structures+W_propulsion+W_equipment_strike+W_useful+W_fuel_strike
W_midfuel_midpayload_a2a = W_structures+W_propulsion+W_equipment_static+W_useful+(0.5*W_payload_a2a)+(0.6*W_fuel_a2a)
W_midfuel_midpayload_strike = W_structures+W_propulsion+W_equipment_static+W_useful+(0.5*W_payload_strike)+(0.6*W_fuel_strike)
# ========================
# xCG Locations
xCG = + origin_offset
xCG_wing = 18.428 + origin_offset
xCG_vt = 28.711 + origin_offset
xCG_fuselage = 18.027 + origin_offset
xCG_ais = 12.194 + origin_offset
xCG_canopy = 4.387 + origin_offset
xCG_en = 16.063 + origin_offset
xCG_fuelsystanks = 20.500 + origin_offset
xCG_furn = 7.026 + origin_offset
xCG_avionics = 10.624 + origin_offset
xCG_pilot = 6.647 + origin_offset
xCG_fuel = 20.500 + origin_offset
xCG_nlg = 7.182 + origin_offset
xCG_mlg = 22.565 + origin_offset
xCG_arrgear = 29.451 + origin_offset
xCG_JDAM = 22.207+ origin_offset
xCG_AIM120C = 23.014 + origin_offset
xCG_AIM9X = 21.421 + origin_offset

xCG_arr = np.array([xCG_wing, xCG_vt, xCG_fuselage, xCG_ais, xCG_canopy, 
                    xCG_en, xCG_fuelsystanks, xCG_furn, xCG_avionics, xCG_pilot,
                      xCG_fuel, xCG_nlg, xCG_mlg, xCG_arrgear, xCG_JDAM, xCG_AIM120C, xCG_AIM9X ])
# print(xCG_arr)
# =========================
# Moments
lump_wing = xCG_wing*(W_wing)
lump_vt = xCG_vt*W_vt
lump_fuselage = xCG_fuselage*(W_fuselage+W_firewall)
lump_furn = xCG_furn*(W_furn+W_flightcont+W_instr)
lump_ais = xCG_ais*(W_ais+W_hydraulics+W_elec+W_acai+W_encontrols+W_encool+W_oilcool+W_starter)
lump_canopy = xCG_canopy*W_canopy
lump_en = xCG_en*(W_en+W_emounts+W_tailpipe+W_ensect)
lump_fuelsystanks = xCG_fuelsystanks*W_fuelsystanks
lump_pylons = (xCG_AIM120C*num_AIM120C*W_py_AIM120C) + (xCG_AIM9X*num_AIM9X*W_py_AIM9X) + (xCG_JDAM*num_JDAM*W_py_JDAM)
lump_nlg = xCG_nlg*(W_nlg+W_catgear+W_handgear)
lump_mlg = xCG_mlg*W_mlg
lump_arrgear = xCG_arrgear*W_arrgear

#
lump_pilot = xCG_pilot*W_pilot
lump_payload_a2a = (xCG_AIM120C*num_AIM120C*W_AIM120C) + (xCG_AIM9X*num_AIM9X*W_AIM9X)
lump_payload_strike = (xCG_JDAM*num_JDAM*W_JDAM) + (xCG_AIM9X*num_AIM9X*W_AIM9X)
lump_fuel_a2a = xCG_fuelsystanks*W_fuel_a2a
lump_fuel_strike = xCG_fuelsystanks*W_fuel_strike

# =================================
lump_mid_payload_a2a = (xCG_AIM120C*0.5*num_AIM120C*W_AIM120C) + (xCG_AIM9X*0.5*num_AIM9X*W_AIM9X)
lump_mid_payload_strike = (xCG_JDAM*0.5*num_JDAM*W_JDAM) + (xCG_AIM9X*0.5*num_AIM9X*W_AIM9X)
# ======================================================================
lump_static = lump_wing+lump_vt+lump_fuselage+lump_furn+lump_ais+lump_canopy+lump_en+lump_fuelsystanks+lump_pylons+lump_nlg+lump_mlg+lump_arrgear
lump_total_a2a = lump_static+lump_payload_a2a+lump_pilot+lump_fuel_a2a
lump_total_strike = lump_static+lump_payload_strike+lump_pilot+lump_fuel_strike
lump_empty_payload_a2a = lump_static+lump_pilot+lump_fuel_a2a
lump_empty_payload_strike = lump_static+lump_pilot+lump_fuel_strike
lump_empty = lump_static +lump_pilot
lump_midmid_a2a = lump_static+lump_pilot+lump_mid_payload_a2a+(0.5*lump_fuel_a2a)
lump_midmid_strike = lump_static+lump_pilot+lump_mid_payload_strike+(0.5*lump_fuel_strike)
print("=====================")

xCG_empty = lump_empty/W_empty
print("xCG empty:", xCG_empty)

xCG_empty_payload_a2a = lump_empty_payload_a2a/W_empty_payload_a2a
print("a2a xCG no payload + fuel:", xCG_empty_payload_a2a)

xCG_empty_payload_strike = lump_empty_payload_strike/W_empty_payload_strike
print("strike xCG no payload + fuel:", xCG_empty_payload_strike)

xCG_total_a2a = lump_total_a2a/W_total_a2a
print("xCG a2a full payload + fuel:", xCG_total_a2a)

xCG_total_strike = lump_total_strike/W_total_strike
print("xCG strike full payload + fuel:", xCG_total_strike)

xCG_mid_a2a = lump_midmid_a2a/W_midfuel_midpayload_a2a
#print("xCG a2a half payload + half fuel:", xCG_mid_a2a)

xCG_mid_strike = lump_midmid_strike/W_midfuel_midpayload_strike
#print("xCG strike half payload + half fuel:", xCG_mid_strike)

LEMAC = 12.844
LEMAC = LEMAC + origin_offset
print("====================")
xCG_empty_percent = (xCG_empty-LEMAC)/MAC * 100
print("xCG_empty_percent:", xCG_empty_percent)

xCG_empty_payload_a2a_percent = (xCG_empty_payload_a2a-LEMAC)/MAC * 100
print("xCG_empty_payload", xCG_empty_payload_a2a_percent)

xCG_empty_payload_strike_percent = (xCG_empty_payload_strike-LEMAC)/MAC * 100
print("xCG_empty_payload", xCG_empty_payload_strike_percent)

xCG_total_a2a_percent = (xCG_total_a2a-LEMAC)/MAC * 100
print("xCG_total_a2a_percent", xCG_total_a2a_percent)

xCG_total_strike_percent = (xCG_total_strike-LEMAC)/MAC * 100
print("xCG_total_strike_percent:", xCG_total_strike_percent)

xCG_strike_empty = xCG_total_strike_percent - xCG_empty_percent
print("strike to empty %:", xCG_strike_empty)

xCG_a2a_empty = xCG_total_a2a_percent - xCG_empty_percent
print("a2a to empty %:", xCG_a2a_empty)

xCG_strike_no_payload = xCG_total_strike_percent - xCG_empty_payload_strike_percent
print("strike to no payload %:", xCG_strike_no_payload)

xCG_a2a_no_payload = xCG_total_a2a_percent - xCG_empty_payload_a2a_percent
print("a2a to no payload %:", xCG_a2a_no_payload)
# Graph

# ==========================
# fuel arrays only
resol = 1000 #resolution (data points)
arr_fuel_a2a = np.linspace(0, W_fuel_a2a, resol)
arr_fuel_strike = np.linspace(0, W_fuel_strike, resol)

arr_fuel_a2a = np.array(arr_fuel_a2a)
arr_fuel_strike = np.array(arr_fuel_strike)

arr_moment_a2a = arr_fuel_a2a*xCG_fuel
arr_moment_strike = arr_fuel_strike*xCG_fuel

# ===========================
# combining fuel arrays into other weights for operations
# ===========================

# Moments
# full payload
M_arr_oper_a2a_full = (lump_static+lump_pilot+lump_payload_a2a)+arr_moment_a2a
M_arr_oper_strike_full = (lump_static+lump_pilot+lump_payload_strike)+arr_moment_strike

# empty payload
M_arr_oper_a2a_empty = (lump_static+lump_pilot)+arr_moment_a2a
M_arr_oper_strike_empty = (lump_static+lump_pilot)+arr_moment_strike

# =========================

# Weights
# full
W_arr_oper_a2a_full = (W_empty+W_pilot+W_payload_a2a)+arr_fuel_a2a
W_arr_oper_strike_full = (W_empty+W_pilot+W_payload_strike)+arr_fuel_strike

# empty payload
W_arr_oper_a2a_empty = (W_empty+W_pilot)+arr_fuel_a2a
W_arr_oper_strike_empty = (W_empty+W_pilot)+arr_fuel_strike

#print(W_arr_oper_a2a_empty[0])
#print(W_arr_oper_strike_empty[0])
#print(W_arr_oper_a2a_empty[-1])
#print(W_arr_oper_strike_empty[-1])

xCG_arr_a2a_full = M_arr_oper_a2a_full/W_arr_oper_a2a_full
xCG_arr_a2a_empty = M_arr_oper_a2a_empty/W_arr_oper_a2a_empty
xCG_arr_strike_full = M_arr_oper_strike_full/W_arr_oper_strike_full
xCG_arr_strike_empty = M_arr_oper_strike_empty/W_arr_oper_strike_empty

plt.title('Change in weight due to fuel burn vs xCG location')
plt.xlabel("xCG (ft)")
plt.ylabel("Weight (lbf)")
#plt.plot(W_arr_oper_a2a_full, xCG_arr_a2a_full, label='A2A full payload')
#plt.plot(W_arr_oper_a2a_empty, xCG_arr_a2a_empty, label='A2A empty payload')
#plt.plot(W_arr_oper_strike_full, xCG_arr_strike_full, label='Strike full payload')
#plt.plot(W_arr_oper_strike_empty, xCG_arr_strike_empty, label='Strike empty payload')
plt.plot(xCG_arr_a2a_full, W_arr_oper_a2a_full, label='A2A full payload')
plt.plot(xCG_arr_a2a_empty, W_arr_oper_a2a_empty, label='A2A empty payload', marker='x')
plt.plot(xCG_arr_strike_full, W_arr_oper_strike_full, label='Strike full payload')
plt.plot(xCG_arr_strike_empty, W_arr_oper_strike_empty, label='Strike empty payload')
plt.legend(loc='best')
plt.grid()
plt.show()