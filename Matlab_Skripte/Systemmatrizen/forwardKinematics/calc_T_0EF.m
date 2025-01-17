function T_0EF = calc_T_0EF(in1,l1,l2)
%calc_T_0EF
%    T_0EF = calc_T_0EF(IN1,L1,L2)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    10-Jan-2025 11:22:53

alpha = in1(1,:);
beta = in1(2,:);
t2 = cos(alpha);
t3 = cos(beta);
t4 = sin(alpha);
t5 = sin(beta);
t6 = t4.*t5;
t7 = t2.*t3;
t8 = t2.*t5;
t9 = t3.*t4;
t10 = -t7;
t11 = t8+t9;
t12 = t6+t10;
T_0EF = reshape([t12,t11,0.0,0.0,-t11,t12,0.0,0.0,0.0,0.0,1.0,0.0,l1.*t4+l2.*t12,l1.*t2+l2.*t11,0.0,1.0],[4,4]);
end
