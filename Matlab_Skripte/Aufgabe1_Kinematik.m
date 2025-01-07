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

%Achse 1:
%DH: a_0 = 0; alpha_0 = 0; d_1 = 0; theta_1 = alpha

%T_12_trans_1 = eye(4);

%T_12_rot_1 = eye(4);

%T_01_trans_2 = eye(4);

T_01 = [cos(alpha), -sin(alpha), 0, 0;
        sin(alpha), cos(alpha), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];                    %=T_01_rot_2


%Achse 2:
%DH: a_1 = l1; alpha_1 = alpha; d_2 = 0; theta_2 = beta
T_12_trans_1 = [1, 0, 0, l1;
              0, 1, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];

T_12_rot_1 = [1, 0, 0, 0;
              0, cos(alpha), -sin(alpha), 0;
              0, sin(alpha), cos(alpha), 0;
              0, 0, 0, 1];


%T_12_trans_2 = eye(4);

T_12_rot_2 = [cos(beta), -sin(beta), 0, 0;
              sin(beta), cos(beta), 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];

T_12 = T_12_trans_1 * T_12_rot_1 * T_12_rot_2;

%Achse 3:
%DH: a_2 = l2; alpha_2 = beta; d_3 = 0; theta_3 = 0
T_23_trans_1 = [1, 0, 0, l2;
              0, 1, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];

T_23_rot_1 = [1, 0, 0, 0;
              0, cos(beta), -sin(beta), 0;
              0, sin(beta), cos(beta), 0;
              0, 0, 0, 1];

%T_23_trans_2 = eye(4);

%T_23_rot_2 = eye(4);

T_23 = T_23_trans_1 * T_23_rot_1;

%Gesamte Transformation:
T_03 = T_01 * T_12 * T_23;  %POSE Endeffektor

disp('T_03 = ');
disp(T_03);
%-------------------------------------------------------------------------------------%
%Ich weiß, dass es falsch ist, habe die DH-Konevtion nicht ganz verstanden, also mache
%nochmal aber richtig ohne DH-Konvention, korrigiere später
%-------------------------------------------------------------------------------------%

T_01 = [cos(alpha), -sin(alpha), 0, 0;
        sin(alpha), cos(alpha), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];

T_12_trans = [1, 0, 0, l1;
              0, 1, 0, 0;
              0, 0, 1, 0;
              0, 0, 0, 1];  %in KOS 1
              
T_12_rot =  [cos(beta), -sin(beta), 0, 0;
             sin(beta), cos(beta), 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];   %rot um z1-Achse

T_12 = T_12_trans * T_12_rot;

T_23 = [1, 0, 0, l2;
        0, 1, 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1];  %in KOS 2 

T_0EF = T_01 * T_12 * T_23;  %POSE Endeffektor

disp('T_0EF = ');
disp(T_0EF);

%Für später ob DH richtig verstanden oder nicht, aktuell noch falsch
delta = T_0EF - T_03;
disp('delta = ');
disp(delta);

%Exportieren von T_0EF als Matlab Function für inverse Kinematik
matlabFunction(T_0EF,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\T_0EF_fcn', 'Vars', {alpha, beta, l1, l2});

%-------------------------------------------------------------------------------------%
%Man kann die symbolic Toolbox nicht in Simulink verwenden...

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
matlabFunction(J,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\J_func', 'Vars', {y, l1, l2});
matlabFunction(J_dot,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\J_dot_func', 'Vars', {y, y_dot, l1, l2});