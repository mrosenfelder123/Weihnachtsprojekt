function M = calc_M(in1,in2)
%calc_M
%    M = calc_M(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    13-Jan-2025 10:11:35

K1_I1 = in2(6,:);
K2_I2 = in2(7,:);
alpha = in1(1,:);
beta = in1(2,:);
l1 = in2(1,:);
l2 = in2(2,:);
m1 = in2(3,:);
m2 = in2(4,:);
t2 = cos(alpha);
t3 = cos(beta);
t4 = sin(alpha);
t5 = sin(beta);
t6 = l1.^2;
t7 = l2.^2;
t8 = alpha.*t2;
t9 = l1.*t2;
t10 = alpha.*t4;
t11 = l1.*t4;
t12 = t4.*t5;
t13 = -t2;
t15 = t2.*t5;
t16 = t3.*t4;
t17 = -t10;
t18 = t4+t8;
t19 = t3.*t13;
t21 = t15+t16;
t20 = t2+t17;
t22 = beta.*t21;
t23 = t12+t19;
t25 = (l2.*t21)./2.0;
t24 = beta.*t23;
t27 = (l2.*t23)./2.0;
t28 = t18+t22;
t29 = t9+t25;
t32 = t22+t23;
t26 = -t24;
t30 = t11+t27;
t31 = t10+t13+t24;
t34 = t21.*t28;
t35 = t23.*t28;
t37 = t21.*t32;
t38 = m2.*t25.*t29;
t42 = t23.*t32;
t33 = t21+t26;
t36 = t21.*t31;
t39 = t23.*t31;
t44 = m2.*t27.*t30;
t40 = -t36;
t41 = t21.*t33;
t43 = t23.*t33;
t46 = t34+t39;
t45 = -t43;
t47 = t35+t40;
t48 = t41+t42;
t49 = t37+t45;
M = reshape([m2.*t29.^2+m2.*t30.^2+K1_I1.*t18.*(t2.*t20+t4.*t18)-K1_I1.*t20.*(t2.*t18-t4.*t20)+(m1.*t2.^2.*t6)./4.0+(m1.*t4.^2.*t6)./4.0+K2_I2.*t28.*t47+K2_I2.*t31.*t46,t38+t44+K2_I2.*t28.*t48+K2_I2.*t31.*t49,t38+t44+K2_I2.*t32.*t47-K2_I2.*t33.*t46,(m2.*t7.*t21.^2)./4.0+(m2.*t7.*t23.^2)./4.0+K2_I2.*t32.*t48-K2_I2.*t33.*t49],[2,2]);
end
