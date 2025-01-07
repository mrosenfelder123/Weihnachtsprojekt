%-------------------------------------------------------------------------------------%
%MSM-Weihnachtsprojekt
%Andreas Paul
%Mat.:3418246
%06.01.2025
%-------------------------------------------------------------------------------------%
clc;
clear;
close all;

%TEILAUFGABE a): 
%Wahl geeigneter Verallgemeinerter Koordinaten

%Verallgemeinerte Koordinaten:
syms alpha beta alpha_dot beta_dot alpha_ddot beta_ddot

y = [alpha; beta];
y_dot = [alpha_dot; beta_dot];
y_ddot = [alpha_ddot; beta_ddot];

%Gravitationsvektor in neg y-Richtung
syms g
I_g = [0; -g; 0];   

%Roboterparameter:
syms l1 l2 m1 m2 I1 I2
%Durchmesser der Gelenke zur bestimmung des wirkenden Reibmoments
syms d1 d2

sys_param = [l1; l2; m1; m2; g; d1; d2; I1; I2];

%Motormomente:	
syms tau1 tau2
tau_m = [tau1; tau2];
sys_input = tau_m;

%TEILAUFGABE b):
%Bestimmung der Bewegungsgleichung

%(i) Bestimmung der Jacobi-Matrix und der kinetischen Energie T

%Transformationsmatritzen aus Aufgabe 1:
%Notiz: I... = Inertialsystem, SP... = Schwerpunkt, EF... = Endeffektor

T_I1 = [cos(alpha), -sin(alpha), 0, 0;
        sin(alpha), cos(alpha), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];

T_12 = [cos(beta), -sin(beta), 0, l1;
        sin(beta),  cos(beta), 0,  0;
        0,          0,         1,  0;
        0,          0,         0,  1];

%Bestimmung der POSE der Schwerpunkte der Körper 1 und 2

%Position des Schwerpunktes des Körpers 1:
T_1SP1 =     [1, 0, 0, l1/2;
              0, 1, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];

T_ISP1 = T_I1 * T_1SP1;

%im Inertialsystem, für E_pot
I_r_SP1 = T_ISP1(1:3, 4);

%im körperfesten Koordinatensystem, für E_kin
K1_r_SP1 = T_1SP1(1:3, 4);

%Orienteirung des Schwerpunktes des Körpers 1:
%Notiz: Berechnung im körperfesten Koordinatensystem, s.d. Trägheitsmoment = konst.
K1_s_ISP1 = [alpha; 0; 0];
I_S_ISP1 = T_ISP1(1:3, 1:3);

%Position des Schwerpunktes des Körpers 2:
T_2SP2 =     [1, 0, 0, l2/2;
              0, 1, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];

T_ISP2 = T_I1 * T_12 * T_2SP2;

%im Inertialsystem, für E_pot
I_r_SP2 = T_ISP2(1:3, 4);

%im körperfesten Koordinatensystem, für E_kin
K2_r_SP2 = T_2SP2(1:3, 4);

%Orienteirung des Schwerpunktes des Körpers 2:
%Notiz: Berechnung im körperfesten Koordinatensystem, s.d. Trägheitsmoment = konst.
K2_s_K1SP2 = [beta; 0; 0];
I_S_ISP2 = T_ISP2(1:3, 1:3);

%Rotationen sind additiv
K2_s_ISP2 = T_12(1:3,1:3).' * K1_s_ISP1 + K2_s_K1SP2;

%Bestimmung der Jacobi-Matrix
I_J_trans_SP1 = jacobian(I_r_SP1, y);
K1_J_rot_SP1 = jacobian(K1_s_ISP1, y);

I_J_trans_SP2 = jacobian(I_r_SP2, y);
K2_J_rot_SP2 = jacobian(K2_s_ISP2, y);

%Bestimmung der Massenmatritzen M1 und M2
%Notiz: I1, I2 bzgl. körperfestem Koordinatensystem!
%FRAGE: Brauche ich Trafo I_S_ISP1)?!?! eig nicht oder?
M1 = m1 * I_J_trans_SP1.' * I_J_trans_SP1 + K1_J_rot_SP1.' * I_S_ISP1 * I1 * I_S_ISP1.' * K1_J_rot_SP1;
M2 = m2 * I_J_trans_SP2.' * I_J_trans_SP2 + K2_J_rot_SP2.' * I_S_ISP2 * I2 * I_S_ISP2.' * K2_J_rot_SP2;

%correct??
M = M1 + M2;

%(ii) Bestimmung der potentiellen Energie U
%via skalarprodukt
U1 = m1 * I_g.' * I_r_SP1;
U2 = m2 * I_g.' * I_r_SP2;

U = U1 + U2;

%(iii) Bestimmung der nicht-konservative Kräfte
%Reibmoment aus viskoser und statischer Reibung

%Notiz: Die sign-Funtion wird approximiert durch (2/pi) * arctan(),
%siehe Reibmoment_TESTFILE.m

%ÜBERARBEITEN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%festgelegte Reibparameter
syms mu_v1 mu_v2 F_s1 F_s2
reib_param = [mu_v1; mu_v2; F_s1; F_s2];


%Gelenk 1
b1 = 20;    %Erhöhung der Steigung von arctan
F_N1 = (m1 + m2) * g;
f_r = mu_c * (2/pi) * atan(b1 * alpha_dot) + f_v * alpha_dot + f_m * exp(-c * abs(alpha_dot)) * sign(alpha_dot);
F_R1 = f_r * F_N1;
tau_R_1 = d1 * F_R1;

%Gelenk 2
b2 = 20;    %Erhöhung der Steigung von arctan
F_N2 = m2 * g;
f_r = mu_c * (2/pi) * atan(b2 * beta_dot) + f_v * beta_dot + f_m * exp(-c * abs(beta_dot)) * sign(beta_dot);
F_R2 = f_r * F_N2;
tau_R_2 = d2 * F_R2;

%Gesamtes Reibmoment
tau_R = [tau_R_1; tau_R_2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(iv) Bestimmung der Bewegungsgleichung in symbolischer Matrixform
%Form: M(y) * y_ddot + D(y, y_dot) * y_dot + G(y) = tau_R + tau_M = q(y, y_dot)

%Berechnung von D(y, y_dot) mit Christoffel-Symbolen
D = sym(zeros(2,2));    %Initialisierung D
for i = 1:2
    for j = 1:2
        for k = 1:2
            D(i,j) = D(i,j) + 0.5 * (diff(M(k,j), y(i)) + diff(M(k,i), y(j)) - diff(M(i,j), y(k))) * y_dot(i);
        end
    end
end

%Brechnung von g(y)
G = jacobian(U, y).';

%Erstellen von Matlab Funktionen für die Teile der Bewegungsgleichungen
matlabFunction(M,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\MassMatrix', 'Vars', {y, sys_param, sys_input});
matlabFunction(D,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\DMatrix', 'Vars', {y, y_dot, sys_param, sys_input, reib_param});
matlabFunction(G,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\GVector', 'Vars', {y, sys_param, sys_input, reib_param});
matlabFunction(tau_R,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\CALCtau_R', 'Vars', {y, y_dot, sys_param, sys_input, reib_param});