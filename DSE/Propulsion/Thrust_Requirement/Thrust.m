clear all;
close all;
clc;

T_req = 105772;
V_0 = 0;
p_0 = 101325;

PI_fan = 1.4;
PI_core = 1;
BPR = 15;

size_arr = 100

mfc_min = 20;
mfc_max = 100;

vfan_min = 20;
vfan_max = 300;

massflow_arr = mfc_min:mfc_max/size_arr:mfc_max;

fan_vel_arr = vfan_min:vfan_max/size_arr:vfan_max;

[mf,vf] = meshgrid(massflow_arr, fan_vel_arr);

vel_core = T_req./mf - vf.*BPR;

surf(mf, vf, vel_core)