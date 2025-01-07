%-------------------------------------------------------------------------------------%
%MSM-Weihnachtsprojekt
%Andreas Paul
%Mat.:3418246
%06.01.2025
%-------------------------------------------------------------------------------------%
clc;
clear;
close all;
%INITIALISIERUNG FÜR SIMULINKMODELL

%System Parameter
%I1, I2 in Modell Function selbst definiert
l1 = 0.16;
l2 = 0.128;
m1 = 0.1817;
m2 = 0.0944;
g  = 9.81;
I1 = 0.0309;
I2 = 0.0045;

%Reib Parameter
%Es können noch b1, b2 eingestellt werden, die die Steigung der 
%atan Funktion, welche genutzt wird zum ersetzten der sign Funktion
%siehe Aufgabe2_Dynamik.m line 133.
mu_v1 =3.843e-6;
mu_v2 = 3.887e-6; 
F_s1  = 8.5e-4;
F_s2  = 3.2e-4;
