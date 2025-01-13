clc;
clear;
close all;

%param def
l1 = 0.16;
l2 = 0.128;
syms alpha alpha_dot alpha_ddot beta beta_dot beta_ddot t
y = [alpha; beta];
y_dot = [alpha_dot; beta_dot];
y_ddot = [alpha_ddot; beta_ddot];

%import POSE of EF and jacobians
T_0EF = calc_T_0EF(y,l1,l2);
J = calc_J(y,l1,l2);
J_dot = calc_J_dot(y,y_dot,l1,l2);

%Position, Geschw, Beschl
r = T_0EF(1:3,4);
r_dot = J * y_dot;
r_ddot = J_dot * y_dot + J * y_ddot;

%get trajektorie
tra = calc_Trajektorie_kartesisch(t);
rd = tra(1:3);
rd_dot = tra(4:6);
rd_ddot = tra(7:9);
%% Berechne yd

%Berechne yd
t_ = 0:0.05:2;
t_steps = length(t_);
yd_alpha = zeros(1,t_steps);
yd_beta = zeros(1,t_steps);

for i = 1:t_steps
rd1 = subs(rd,t,t_(i));
eq1 = r == rd1;

yd1 = solve(eq1, alpha, beta);
alpha_sol = double(yd1.alpha(1));
beta_sol = double(yd1.beta(1));
yd_alpha(i) = alpha_sol;
yd_beta(i) = beta_sol;
end

%fitting Polynom 3. Ordnung zu yd
%-->Korrektur fitting Ordnung 5, da 3 nicht gut genug
p_alpha = polyfit(t_, yd_alpha, 5);

yd_alpha_sol = p_alpha(1) * t^5 + p_alpha(2) * t^4 + p_alpha(3) * t^3 + p_alpha(4) * t^2 +...
    p_alpha(5) * t + p_alpha(6);

p_beta = polyfit(t_, yd_beta, 5);

yd_beta_sol = p_beta(1) * t^5 + p_beta(2) * t^4 + p_beta(3) * t^3 + p_beta(4) * t^2 +...
    p_beta(5) * t + p_beta(6);

yd = [yd_alpha_sol; yd_beta_sol];

%% Berechne yd_dot

yd_dot_alpha_dot = zeros(1,t_steps);
yd_dot_beta_dot = zeros(1,t_steps);

for i = 1:t_steps
rd_dot1 = subs(rd_dot,t,t_(i));

%subs yd in r_dot
alphad1 = subs(yd_alpha_sol,t,t_(i));
betad1 = subs(yd_beta_sol,t,t_(i));
yd2 = [alphad1; betad1];
r_dot1 = subs(r_dot, y,yd2);

eq2 = r_dot1 == rd_dot1;
yd_dot1 = solve(eq2, alpha_dot, beta_dot);
alpha_dot_sol = double(yd_dot1.alpha_dot(1));
beta_dot_sol = double(yd_dot1.beta_dot(1));
yd_dot_alpha_dot(i) = alpha_dot_sol;
yd_dot_beta_dot(i) = beta_dot_sol;
end

%fitting Ordnung 5, da 3 nicht gut genug
p_alpha_dot = polyfit(t_, yd_dot_alpha_dot, 5);

yd_dot_alpha_dot_sol = p_alpha_dot(1) * t^5 + p_alpha_dot(2) * t^4 + p_alpha_dot(3) * t^3 + p_alpha_dot(4) * t^2 +...
    p_alpha_dot(5) * t + p_alpha_dot(6);

p_beta_dot = polyfit(t_, yd_dot_beta_dot, 5);

yd_dot_beta_dot_sol = p_beta_dot(1) * t^5 + p_beta_dot(2) * t^4 + p_beta_dot(3) * t^3 + p_beta_dot(4) * t^2 +...
    p_beta_dot(5) * t + p_beta_dot(6);

yd_dot = [yd_dot_alpha_dot_sol; yd_dot_beta_dot_sol];

%% Berechne yd_ddot

yd_ddot_alpha_ddot = zeros(1,t_steps);
yd_ddot_beta_ddot = zeros(1,t_steps);

for i = 1:t_steps
rd_ddot1 = subs(rd_ddot,t,t_(i));

%subs yd yd_dot in r_ddot
alphad2 = subs(yd_alpha_sol,t,t_(i));
betad2 = subs(yd_beta_sol,t,t_(i));
alpha_dotd2 = subs(alpha_dot_sol,t,t_(i));
betad_dot2 = subs(beta_dot_sol,t,t_(i));

yd3 = [alphad2; betad2];
yd_dot3 = [alpha_dotd2; betad_dot2];

r_ddot1 = subs(r_ddot, y,yd3);
r_ddot1 = subs(r_ddot1, y_dot,yd_dot3);

eq3 = r_ddot1 == rd_ddot1;
yd_ddot1 = solve(eq3, alpha_ddot, beta_ddot);
alpha_ddot_sol = double(yd_ddot1.alpha_ddot(1));
beta_ddot_sol = double(yd_ddot1.beta_ddot(1));
yd_ddot_alpha_ddot(i) = alpha_ddot_sol;
yd_ddot_beta_ddot(i) = beta_ddot_sol;
end

%fitting Ordnung 5, da 3 nicht gut genug
p_alpha_ddot = polyfit(t_, yd_ddot_alpha_ddot, 5);

yd_ddot_alpha_ddot_sol = p_alpha_ddot(1) * t^5 + p_alpha_ddot(2) * t^4 + p_alpha_ddot(3) * t^3 + p_alpha_ddot(4) * t^2 +...
    p_alpha_ddot(5) * t + p_alpha_ddot(6);

p_beta_ddot = polyfit(t_, yd_ddot_beta_ddot, 5);

yd_ddot_beta_ddot_sol = p_beta_ddot(1) * t^5 + p_beta_ddot(2) * t^4 + p_beta_ddot(3) * t^3 + p_beta_ddot(4) * t^2 +...
    p_beta_ddot(5) * t + p_beta_ddot(6);

yd_ddot = [yd_ddot_alpha_ddot_sol; yd_ddot_beta_ddot_sol];

%% Plot die Funktionen zur Überprüfung

% Evaluate yd_alpha_sol at the given time points
yd_alpha_sol_points = polyval(p_alpha, t_);

% Plot yd_alpha_sol against some points of yd_alpha
figure;
plot(t_, yd_alpha, 'bo', 'DisplayName', 'Original Data');
hold on;
plot(t_, yd_alpha_sol_points, 'r-', 'DisplayName', 'Fitted Curve', 'LineWidth', 1);
xlabel('Time $t$ (s)', 'Interpreter', 'latex');
ylabel('$\alpha$ (rad)', 'Interpreter', 'latex');
title('Fitting Curve to $\alpha$ over Time', 'Interpreter', 'latex');
legend('show');
grid on;

% Evaluate yd_dot_alpha_dot_sol at the given time points
yd_dot_alpha_dot_sol_points = polyval(p_alpha_dot, t_);

% Plot yd_dot_alpha_dot_sol against some points of yd_dot_alpha_dot
figure;
plot(t_, yd_dot_alpha_dot, 'bo', 'DisplayName', 'Original Derivative Data');
hold on;
plot(t_, yd_dot_alpha_dot_sol_points, 'r-', 'DisplayName', 'Fitted Derivative Curve', 'LineWidth', 1);
xlabel('Time $t$ (s)', 'Interpreter', 'latex');
ylabel('$\dot{\alpha}$ (rad/s)', 'Interpreter', 'latex');
title('Fitting Derivative Curve to $\dot{\alpha}$ over Time', 'Interpreter', 'latex');
legend('show');
grid on;

% Evaluate yd_ddot_alpha_ddot_sol at the given time points
yd_ddot_alpha_ddot_sol_points = polyval(p_alpha_ddot, t_);

% Plot yd_ddot_alpha_ddot_sol against some points of yd_ddot_alpha_ddot
figure;
plot(t_, yd_ddot_alpha_ddot, 'bo', 'DisplayName', 'Original Second Derivative Data');
hold on;
plot(t_, yd_ddot_alpha_ddot_sol_points, 'r-', 'DisplayName', 'Fitted Second Derivative Curve', 'LineWidth', 1);
xlabel('Time $t$ (s)', 'Interpreter', 'latex');
ylabel('$\ddot{\alpha}$ (rad/s^2)', 'Interpreter', 'latex');
title('Fitting Second Derivative Curve to $\ddot{\alpha}$ over Time', 'Interpreter', 'latex');
legend('show');
grid on;

%% EXPORT ALS FUNKTIONEN
Trajectory = [yd; yd_dot; yd_ddot];
disp('Trajektorie = ')
disp(Trajectory);
matlabFunction(Trajectory,'File', ...
    'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\Trajektorien\calc_Trajektorie_verallg','Vars', {t});