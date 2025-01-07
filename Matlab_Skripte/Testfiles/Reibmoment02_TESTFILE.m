% Parameter festlegen
b1 = 20;    % Erhöhung der Steigung von arctan
F_s1 = 1;   % Beispielwert für F_s1
mu_v1 = 0.1; % Beispielwert für mu_v1

% Wertebereich für alpha_dot
alpha_dot = linspace(-10, 10, 1000);

% Funktion tau_R_1 berechnen
tau_R_1 = (2/pi) * atan(b1 * alpha_dot) .* (F_s1 + abs(alpha_dot) * mu_v1);

% Plot erstellen
figure;
plot(alpha_dot, tau_R_1);
xlabel('\alpha_{dot}');
ylabel('\tau_{R1}');
title('Reibmoment \tau_{R1} über \alpha_{dot}');
grid on;