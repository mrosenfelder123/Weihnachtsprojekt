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
T_0EF = T_0EF_fcn(alpha,beta,l1,l2);
J = Jfunc(y,l1,l2);
J_dot = J_dot_func(y,y_dot,l1,l2);

%Position, Geschw, Beschl
r = T_0EF(1:3,4);
r_dot = J * y_dot;
r_ddot = J_dot * y_dot + J * y_ddot;

%get trajektorie
tra = CALCtrajektorie(t);
rd = tra(1:3);
rd_dot = tra(4:6);
rd_ddot = tra(7:9);

%Berechne yd
t_ = 0:0.1:2;
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
p_alpha = polyfit(t_, yd_alpha, 3);

yd_alpha_sol = p_alpha(1) * t^3 + p_alpha(2) * t^2 + p_alpha(3) * t + p_alpha(4);

p_beta = polyfit(t_, yd_beta, 3);

yd_beta_sol = p_beta(1) * t^3 + p_beta(2) * t^2 + p_beta(3) * t + p_beta(4);

yd = [yd_alpha_sol; yd_beta_sol];


%plot um zu überprüfen ob Ergebnis stetig
yd_alpha_fitted = polyval(p_alpha, t_);
figure;
plot(t_, yd_alpha, 'bo', 'DisplayName', 'Original Data');
hold on;
plot(t_, yd_alpha_fitted, 'r-', 'DisplayName', 'Fitted Curve', 'LineWidth', 1);
xlabel('Time $t$ (s)', 'Interpreter', 'latex');
ylabel('$\alpha$ (rad)', 'Interpreter', 'latex');
title('Fitting Curve to $\alpha$ over Time', 'Interpreter', 'latex');
legend('show');
grid on;

%Berechne yd_dot


%HIER WEITER

%ACHTE DRAUF, dass hier J und rd_dot subs im for loop drinne steht, sonst
%kommt scheiße

%ALTER FOR LOOP von oben, sprich anpassen
% for i = 1:t_steps
% rd1 = subs(rd,t,t_(i));
% eq1 = r == rd1;
% 
% yd1 = solve(eq1, alpha, beta);
% alpha_sol = double(yd1.alpha(1));
% beta_sol = double(yd1.beta(1));
% yd_alpha(i) = alpha_sol;
% yd_beta(i) = beta_sol;
% end


%EXPORT ALS FUNKTIONEN

%anpassen
% matlabFunction(M,'File', 'D:\MASTER\Semester3\MSM\Weihnachtsprojekt\Matlab_Skripte\Systemmatrizen\MassMatrix', 'Vars', {y, sys_param, sys_input});