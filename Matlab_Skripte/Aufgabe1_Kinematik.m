%-------------------------------------------------------------------------------------%
%MSM-Weihnachtsprojekt
%Andreas Paul
%Mat.:3418246
%06.01.2025
%-------------------------------------------------------------------------------------%
clc;
clear;
close all;
%Transformationsmatrizen für die Kinematik des Roboters

%Verallgemeinerte Koordinaten:
syms alpha beta

%für Jacobians
syms alpha_dot beta_dot
y = [alpha; beta];
y_dot = [alpha_dot; beta_dot];

%Roboterparameter:
syms l1 l2

%Berechung der Transformationsmatritzen
T_01 = [sin(alpha), -cos(alpha), 0, 0;
        cos(alpha),  sin(alpha), 0, 0;
        0,           0,          1, 0;
        0,           0,          0, 1];

matlabFunction(T_01,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\forwardKinematics\calc_T_01', ...
               'Vars', {y,l1,l2});

T_12_trans = [1, 0, 0, l1;
              0, 1, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];  %in KOS 1
              
T_12_rot =  [sin(beta), -cos(beta), 0, 0;
             cos(beta),  sin(beta), 0, 0;
             0,          0,         1, 0;
             0,          0,         0, 1];   %rot um z1-Achse

T_12 = T_12_trans * T_12_rot;

T_02 = T_01 * T_12;

matlabFunction(T_02,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\forwardKinematics\calc_T_02', ...
               'Vars', {y,l1,l2});

T_23 = [1, 0, 0, l2;
        0, 1, 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];  %=T_2EF; in KOS 2 

T_0EF = T_01 * T_12 * T_23;  %POSE Endeffektor

disp('T_0EF = ');
disp(T_0EF);

matlabFunction(T_0EF,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\forwardKinematics\calc_T_0EF', ...
               'Vars', {y,l1,l2});

%-------------------------------------------------------------------------------------%
%Man kann die symbolic Toolbox nicht in Simulink verwenden...
%-------------------------------------------------------------------------------------%

%exportiere J, J_dot --> Gedankengang siehe InverseKinematics_TESTFILE.m
r = T_0EF(1:3,4);

%Geschwindigkeit
%Notiz: r ist in inertial coordiantes? KA deutsch lol, deswegen physical
%derr. gleich numerischer derr. 
%und nutze Kettenregel: r_dot = (dr/dy)*y_dot
J = jacobian(r,y);

%Beschleunigung
%erneut Kettenregel:
%r_ddot = J_dot*y_dot + J*y_ddot
%Da J eine Matrix ist, kann man nicht einfach mit jacobian Befehl arbeiten,
%muss Elementweise erledigt werden.
J_dot = sym(zeros(size(J))); 
for i = 1:size(J, 1)
    for j = 1:size(J, 2)
        J_dot(i, j) = jacobian(J(i, j), y) * y_dot;
    end
end

%Exportieren als Funktionen
matlabFunction(J,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\forwardKinematics\calc_J', ...
               'Vars', {y, l1, l2});
matlabFunction(J_dot,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\forwardKinematics\calc_J_dot', ...
               'Vars', {y, y_dot, l1, l2});