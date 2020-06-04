

#%% ---------------------- Constants ----------------------
p_amb_to = 101325 #ambient pressure during take-off
p_amb_cr = 22632.1 #ambient pressure during cruise
T_amb_to = 288.15 #ambient temperature during take-off
T_amb_cr = 216.65 #ambient temperature during cruise
v_0_to = 0 #freestream velocity during take-off
v_0_cr = 230 #freestream velocity during cruise

cp_air = 1000 #cp value during cold conditions
cp_gas = 1150 #cp value during hot conditions
k_air = 1.4 #k value during cold conditions
k_gas = 1.333 #k value during hot conditions

#%% ---------------------- Inputs ----------------------
BPR = 15 #bypass ratio
M_core_to = 53 #massflow through the core during take-off
M_core_cr = 20 #massflow through the core during cruise
T_comb_to = 1250 #turbine inlet temperature during take-off
T_comb_cr = 1200 #turbine inlet temperature during cruise

PI_fan = 1.1 #pressure ratio fan
PI_HPC  = 20 #pressure ratio high pressure compressor
PI_LPC = 2 #pressure ratio low pressure compressor

#expected values
Thrust_to = 105773 #required thrust during take-off
Thrust_cr = 12859 #required thrust during cruise
V_fan_e_to = 116 #required fan exhaust velocity during take-off
V_fan_e_cr = 263 #required fan exhaust velocity during cruise
V_core_e_to = 255.72 #required core exhaust velocity during take-off
V_core_e_cr = 377.96 #required core exhaust velocity during cruise
