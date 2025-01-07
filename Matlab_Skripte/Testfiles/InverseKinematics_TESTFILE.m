rd = [0; -0.228; 0];
rd_dot = [0; -0.228; 0];
rd_ddot = [0; -0.228; 0];

%Initialisierung
syms alpha beta alpha_dot beta_dot alpha_ddot beta_ddot
y = [alpha; beta];
y_dot = [alpha_dot; beta_dot];
y_ddot = [alpha_ddot; beta_ddot];
T_0EF = T_0EF_fcn(alpha, beta, l1, l2);


r = T_0EF(1:3,4);


%Ort
eq1 = r == rd;

q = solve(eq1, alpha, beta);
alpha_sol = double(q.alpha(1));
beta_sol = double(q.beta(1));
y_new = [alpha_sol; beta_sol];

%Geschwindigkeit
%Notiz: r ist in inertial coordiantes? KA deutsch lol, deswegen physical
%derr. gleich numerischer derr. 
%und nutze Kettenregel: r_dot = (dr/dy)*y_dot
J = jacobian(r,y);
r_dot = J * y_dot;
r_dot = subs(r_dot, y, y_new);

eq2 = r_dot == rd_dot;

q_dot = solve(eq2, alpha_dot, beta_dot);
alpha_dot_sol = double(q_dot.alpha_dot(1));
beta_dot_sol = double(q_dot.beta_dot(1));
y_dot_new = [alpha_dot_sol; beta_dot_sol];

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

r_ddot = J_dot * y_dot + J * y_ddot;
r_ddot = subs(r_ddot, [y;y_dot], [y_new; y_dot_new]);

eq3 = r_ddot == rd_ddot;

q_ddot = solve(eq3, alpha_ddot, beta_ddot);
alpha_ddot_sol = double(q_ddot.alpha_ddot(1));
beta_ddot_sol = double(q_ddot.beta_ddot(1));
y_ddot_new = [alpha_ddot_sol; beta_ddot_sol];

yd = [y_new; y_dot_new; y_ddot_new];
