%-------------------------------------------------------------------------------------%
%MSM-Weihnachtsprojekt
%Andreas Paul
%Mat.:3418246
%06.01.2025
%-------------------------------------------------------------------------------------%
clear;
close all;
clc;

%vom Endeffektor anzusteuernde Positionen aus Aufgabenstellung
I_pEF_0 = [0; -0.228; 0];
I_pEF_d = [4*sqrt(6) - 4*sqrt(2) - 10; -4*sqrt(6) - 4*sqrt(2) - 10*sqrt(3)   ; 0]*(1/125); %in m

%Es soll eine ruckbegrenze Trajektorie gewählt werden, daher muss Beschleunigung stetig sein
%-->Wähle cubic-function
t=2;    %Anfahrtszeit

p_x = cubic_func(0,I_pEF_0(1),0,t,I_pEF_d(1),0);
p_y = cubic_func(0,I_pEF_0(2),0,t,I_pEF_d(2),0);

syms t
trajektorie = [p_x(1) + p_x(2)*t + p_x(3)*t.^2 + p_x(4)*t.^3;
               p_y(1) + p_y(2)*t + p_y(3)*t.^2 + p_y(4)*t.^3;
               0];
matlabFunction(trajektorie,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\M', 'Vars', {t});



%HILFSFUNTIONEN
%CODE AUS ILIAS, Begleitmaterialien

% ** Cubic Polynomial Traject0ries **
function [a] = cubic_func(t0,q0,v0,tf,qf,vf)
% Syntax
% ======
%   a=cubic_func(t0,q0,v0,tf,qf,vf)
% Description
% ===========
%   Computes the polynomial coefficients of a cubic polynominal
%   and plots the polynominal, its derivate and second derivative
%
% Input Arguments
% ===============
%
% -t0:   'start time'
% -q0:   'postion at start time'
% -v0:   'velocity at start time'
% -tf:   'final time'
% -qf:   'postion at final time'
% -vf:   'velocity at final time'
%
% Output Arguments
% ================
%   -a: coefficients of polynominal
%--------------------------------------------------------------------------
% Vorlesung Modellierung und Simulation in der Mechatronik
% Copyright ITM University of Stuttgart, <a href="matlab:
%  web('http://www.itm.uni-stuttgart.de')">www.itm.uni-stuttgart.de</a>
%--------------------------------------------------------------------------

% Coefficient matrix for cubic traject0ry and its derivative
% at initial and final joint values.
A = [1,  t0,  t0^2, t0^3;   ...
    0,  1,   2*t0, 3*t0^2; ...
    1,  tf,  tf^2, tf^3;   ...
    0,  1,   2*tf, 3*tf^2];
% Vect0r of intial and final joint positions and velocities
b = [q0; v0; qf; vf];
% Compute coefficients of traject0ry polynomial using
% notion of a = inv(A)*b, but using Gaussian Elimination
a = A\b;
% Evaluate cubic polynomial at  500 times
t = t0:(tf-t0)/500:tf;
q     = a(1) + a(2)*t + a(3)*t.^2 + a(4)*t.^3;
qdot  = a(2) + 2*a(3)*t + 3*a(4)*t.^2;
qddot = 2*a(3) + 6*a(4)*t;

% Plot the Trajectories
figure(1);

subplot(1,3,1);
plot(t,q,'r-.','LineWidth',2);

leg=get(legend,'String');
leg{end+1}='q x^3';
legend(leg);


xlabel('time (sec)'); ylabel('Joint Angle [rad]');
% title('Traject0ry using Cubic Polynomial');
hold on;

subplot(1,3,2);
plot(t,qdot,'r-.','LineWidth',2);


leg=get(legend,'String');
leg{end+1}='dq/dt x^3';
legend(leg);
xlabel('time (sec)'); ylabel('Veclocity [rad/s]');
% title('Traject0ry using Cubic Polynomial');
hold on;

subplot(1,3,3);
plot(t,qddot,'r-.','LineWidth',2);

leg=get(legend,'String');
leg{end+1}='d^2q/dt^2 x^3';
legend(leg);
xlabel('time (sec)'); ylabel('Acceleration [rad/s^2]');
% title('Traject0ry using Cubic Polynomial');
hold on;
end


%MIGHT NOT NECC
function[]=LSPB_func(t0,q0,tf,qf,V)
% Syntax
% ======
%   a=LSPB_func(t0,q0,tf,qf,V)
% Description
% ===========
%   Computes the polynomial coefficients of a cubic polynominal
%   and plots the polynominal, its derivate and second derivative
%
% Input Arguments
% ===============
%
% -t0:   'start time'
% -q0:   'postion at start time'
% -tf:   'final time'
% -qf:   'postion at final time'
% % -tb:   'blend time'
% -V:    'constant velocity during Blend time
%
% Output Arguments
% ================
%   -a: coefficients of polynominal
%--------------------------------------------------------------------------
% Vorlesung Modellierung und Simulation in der Mechatronik
% Copyright ITM University of Stuttgart, <a href="matlab:
%  web('http://www.itm.uni-stuttgart.de')">www.itm.uni-stuttgart.de</a>
%--------------------------------------------------------------------------
tb = (q0 - qf + V*tf)/V;
% check that V is within limits
Vmin = (qf - q0)/tf;
% if (V < Vmin || V > 2*Vmin)
%     display(['V = ',num2str(V), ' is not within limits',...
%         '(',num2str(Vmin),', ',num2str(2*Vmin),')']);
%     display('LSPB will not be correct!');
% end;
alpha=V/tb;
a(1) = q0; a(2) = 0; a(3) = alpha/2;


b(1) = qf - alpha*tf^2/(2); b(2) = alpha*tf; b(3) = -alpha/2;

% ** Linear Segments with Parabolic Blends (PSPB) **
% Begin with unshifted version on t = [0, 7].
t = t0:(tf-t0)/500:tf;

q = (a(1) + a(3)*t.^2).*(t<=tb) + ...
    ((qf + q0 - V*tf)/2 + V*t).*((t>tb)-(t>=(tf-tb))) + ...
    (b(1) + b(2)*t + b(3)*t.^2).*(t>(tf-tb));
qdot = (a(2) + 2*a(3)*t).*(t<=tb) + ...
    V.*((t>tb)-(t>=(tf-tb))) + ...
    (b(2) + 2*b(3)*t).*(t>(tf-tb));
qddot = 2*a(3)*(t<=tb) + ...
    0*((t>tb)-(t>=(tf-tb))) + ...
    2*b(3)*(t>(tf-tb));


figure(1);
subplot(1,3,1);
plot(t,q,'m:','LineWidth',2);
leg=get(legend,'String');
leg{end+1}='q LSPB';
%subplot(2,2,1);

xlabel('time (sec)'); ylabel('Joint Angle [rad]');
% title('Po');
hold on;

subplot(1,3,2);
plot(t,qdot,'m:','LineWidth',2);
leg=get(legend,'String');
leg{end+1}='dq/dt LSPB';
legend(leg);
xlabel('time (sec)'); ylabel('Veclocity [rad/s]');
hold on;
% title('Trajectory using LSPB');
hold;
subplot(1,3,3);
plot(t,qddot,'m:','LineWidth',2);
leg=get(legend,'String');
leg{end+1}='d^2q/dt^2 LSPB';legend(leg);
xlabel('time (sec)'); ylabel('Acceleration [rad/s^2]');
hold on;
end