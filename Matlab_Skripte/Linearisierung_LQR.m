%-------------------------------------------------------------------------------------%
%MSM-Weihnachtsprojekt
%Andreas Paul
%Mat.:3418246
%06.01.2025
%-------------------------------------------------------------------------------------%
clc;
clear;
close all;


%Selbe Def wie in Aufgabe 2
syms alpha beta alpha_dot beta_dot alpha_ddot beta_ddot

y = [alpha; beta];
y_dot = [alpha_dot; beta_dot];

%state
x = [y; y_dot];

l1 = 0.16;
l2 = 0.128;
m1 = 0.1817;
m2 = 0.0944;
g  = 9.81;
I1 = 0.0309;
I2 = 0.0045;

mu_v1 = 3.843e-6;
mu_v2 = 3.887e-6; 
F_s1  = 8.5e-4;
F_s2  = 3.2e-4;

%parameter Definiton
%sys_param = [l1; l2; m1; m2; g; K1_I1; K2_I2];
sys_param = [l1; l2; m1; m2; g; I1; I2];
%reib_param = [mu_v1; mu_v2; F_s1; F_s2];
reib_param = [mu_v1; mu_v2; F_s1; F_s2];

M = calc_M(y, sys_param);	