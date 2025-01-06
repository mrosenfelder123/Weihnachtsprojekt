clear;
close all;
clc;
% Parameter festlegen
mu_c = 0.5;
f_v = 0.1;
f_m = 0.05;
c = 1.0;

% Parameter für die Steilheit der arctan-Kurve
b = 20; % Erhöhter Verschiebungsfaktor für steileren Übergang

% Wertebereich für beta_dot
beta_dot = linspace(-10, 10, 1000);

% Reibmoment berechnen
f_r = mu_c * sign(beta_dot) + f_v * beta_dot + f_m .* exp(-c * abs(beta_dot)) .* sign(beta_dot);

% Approximation der steilen Stelle mit arctan hinzufügen
%sign(beta_dot) = 2/pi * atan(b * beta_dot), 2/pi, sodass der Wertebereich von -1 bis 1 geht
f_r2 = mu_c * (2/pi) * atan(b * beta_dot) + f_v * beta_dot + f_m .* exp(-c * abs(beta_dot)) .* sign(beta_dot);

% Plot erstellen
figure;
plot(beta_dot, f_r, 'b', 'DisplayName', 'Original f_r');
hold on;
plot(beta_dot, f_r2, 'r--', 'DisplayName', 'Approximation mit arctan');
xlabel('\beta_{dot}');
ylabel('f_r');
title('Reibmoment f_r über \beta_{dot}');
legend;
grid on;