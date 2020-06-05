from Core_equations import*

#%% ---------------------- Constants ----------------------
p_amb_to = 101325 #ambient pressure during take-off
p_amb_cr = 22632.1 #ambient pressure during cruise
T_amb_to = 288.15 #ambient temperature during take-off
T_amb_cr = 216.65 #ambient temperature during cruise
v_0_to = 0 #freestream velocity during take-off
v_0_cr = 230.1382 #freestream velocity during cruise

Cp_air = 1000 #cp value during cold conditions
Cp_gas = 1150 #cp value during hot conditions
k_air = 1.4 #k value during cold conditions
k_gas = 1.333 #k value during hot conditions
LHV_f = 4.31E7 #lower heating value fuel

#%% ---------------------- Inputs ----------------------
BP_ratio = 15 #bypass ratio
M_core_to = 53 #massflow through the core during take-off
M_core_cr = 15.9 #massflow through the core during cruise
TIT_to = 1250 #turbine inlet temperature during take-off
TIT_cr = 1350 #turbine inlet temperature during cruise
P_bli = 1.64E6

p_ratio_fan = 1.2 #pressure ratio fan
p_ratio_HPC  = 21 #pressure ratio high pressure compressor
p_ratio_LPC = 2.76 #pressure ratio low pressure compressor
p_ratio_comb = 0.99 #pressure ratio combustion chamber

efficiency_inlet = 0.97
efficiency_fan = 0.91
efficiency_LPC = 0.94
efficiency_HPC = 0.92
efficiency_HPT = 0.94
efficiency_LPT = 0.95
efficiency_comb = 0.99
efficiency_exhaust = 0.95

#expected values
Thrust_to = 105773 #required thrust during take-off
Thrust_cr = 12859 #required thrust during cruise
V_fan_e_to = 116 #required fan exhaust velocity during take-off
V_fan_e_cr = 263 #required fan exhaust velocity during cruise
V_core_e_to = 255.72 #required core exhaust velocity during take-off
V_core_e_cr = 377.96 #required core exhaust velocity during cruise

#%% ---------------------- calculation ----------------------

def Engine(BPR,M_core,TIT,PI_fan,PI_HPC,PI_LPC,PI_comb,eff_inlet,eff_fan,eff_LPC,eff_HPC,eff_HPT,eff_LPT,eff_comb,eff_exh,p_amb,T_amb,v_0,cp_air,cp_gas,K_air,K_gas,Power_BLI,LHV):
    """
    Engine thrust calculation
    Outputs the total thrust of the engine
    Inputs: bypass ratio, massflow core, temp inlet turbine, pressure ratio fan, pressure ratio HPC,
            pressure ratio LPC, pressure ratio combustion, efficiency inlet, efficiency fan,
            efficiency LPC, efficiency HPC, efficiency HPT, efficiency LPT, efficency combustion,
            efficiency exhaust, ambient pressure, ambient temperature, inlet velocity, cp_air, cp_gas,
            K_air, K_gas, Power required for BLI, Lower heating value fuel
    """
    M_fan = M_core*BPR
    T_t2, p_t2 = inlet(T_amb,p_amb,eff_inlet,v_0,cp_air,K_air) #inlet
    T_t21, p_t21, Power_fan = compressor(T_t2,p_t2,eff_fan,PI_fan,M_fan,cp_air,K_air) #fan
    T_t24, p_t24, Power_LPC = compressor(T_t21,p_t21,eff_LPC,PI_LPC,M_core,cp_air,K_air) #LPC
    T_t3, p_t3, Power_HPC = compressor(T_t24,p_t24,eff_HPC,PI_HPC,M_core,cp_air,K_air) #HPC

    Power_HPT = Power_HPC #power provided by HPT
    Power_LPT = Power_LPC + Power_fan + Power_BLI #power provided by LPT

    mf_f, p_t4, M_core = combustion(T_t3,TIT,M_core,LHV,cp_gas,eff_comb,p_t3,PI_comb) #combustion chamber
    T_t45,p_t45, PI_HPT = turbine(Power_HPT,M_core,cp_gas,TIT,p_t4,eff_HPT,K_gas) #HPT
    T_t5,p_t5, PI_LPT = turbine(Power_LPT,M_core,cp_gas,T_t45,p_t45,eff_LPT,K_gas) #LPT

    T_e_fan, v_e_fan = exhaust(p_amb,T_t21,p_t21,eff_exh,K_air,cp_air) #exhaust fan
    T_e_core, v_e_core = exhaust(p_amb,T_t5,p_t5,eff_exh,K_gas,cp_gas) #exhaust core

    thrust_core = thrust_exhaust(v_e_core,v_0,M_core) #thrust core
    thrust_fan = thrust_exhaust(v_e_fan,v_0,M_fan) #thrust fan
    thrust = thrust_core + thrust_fan #total thrust
    return thrust

print(Engine(BP_ratio,M_core_cr,TIT_cr,p_ratio_fan,p_ratio_HPC,p_ratio_LPC,p_ratio_comb,efficiency_inlet,efficiency_fan,efficiency_LPC,efficiency_HPC,efficiency_HPT,efficiency_LPT,efficiency_comb,efficiency_exhaust,p_amb_cr,T_amb_cr,v_0_cr,Cp_air,Cp_gas,k_air,k_gas,P_bli,LHV_f))