function D = D(in1,in2,in3,in4,in5)
%D
%    D = D(IN1,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 24.1.
%    06-Jan-2025 22:32:24

I2 = in3(6,:);
alpha = in1(1,:);
alpha_dot = in2(1,:);
beta = in1(2,:);
beta_dot = in2(2,:);
l1 = in3(1,:);
l2 = in3(2,:);
m2 = in3(4,:);
t2 = cos(alpha);
t3 = cos(beta);
t4 = sin(alpha);
t5 = sin(beta);
t6 = l1.*t2;
t7 = alpha.*t5;
t8 = l1.*t4;
t9 = t4.*t5;
t10 = t2.*t3;
t11 = t2.*t5;
t12 = t3.*t4;
t14 = -t10;
t15 = t7-1.0;
t16 = t11+t12;
t17 = t9+t14;
t18 = t3.*t16;
t19 = t5.*t16;
t23 = t7.*t16;
t24 = (l2.*t16)./2.0;
t39 = t15.*t16;
t20 = t18.*2.0;
t21 = t19.*2.0;
t22 = alpha.*t18;
t25 = t3.*t17;
t26 = -t18;
t27 = t5.*t17;
t34 = t7.*t17;
t35 = (l2.*t17)./2.0;
t40 = t8+t24;
t41 = t15.*t17;
t29 = t25.*2.0;
t30 = t27.*2.0;
t31 = alpha.*t25;
t33 = t22.*-2.0;
t36 = -t27;
t42 = -t41;
t44 = t19+t25;
t50 = m2.*t35.*t40;
t53 = m2.*t24.*(t6-t35);
t37 = -t30;
t38 = alpha.*t29;
t45 = t18+t36;
t46 = t21+t29;
t48 = t31+t39;
t54 = I2.*t16.*t44;
t56 = I2.*t17.*t44;
t71 = -I2.*t17.*(t22+t42);
t72 = I2.*t16.*(t22+t42);
t74 = t33+t34+t41;
t75 = I2.*t17.*(t22+t42);
t79 = t22+t42+t44;
t47 = t20+t37;
t55 = I2.*t16.*t45;
t57 = I2.*t16.*t46;
t59 = I2.*t17.*t45;
t60 = I2.*t17.*t46;
t66 = I2.*t16.*t48;
t68 = I2.*t17.*t48;
t70 = t23+t38+t39;
t78 = I2.*t16.*t74;
t80 = t26+t27+t48;
t82 = I2.*t17.*t74;
t83 = I2.*t16.*t79;
t85 = I2.*t17.*t79;
t58 = I2.*t16.*t47;
t61 = -t57;
t62 = I2.*t17.*t47;
t64 = -t59;
t65 = -t60;
t69 = -t66;
t76 = I2.*t16.*t70;
t77 = I2.*t17.*t70;
t84 = I2.*t16.*t80;
t86 = I2.*t17.*t80;
t88 = t55+t56;
t92 = t66+t71;
t96 = (t3.*(t68+t72))./2.0;
t63 = -t58;
t81 = -t77;
t87 = -t84;
t89 = t54+t64;
t90 = (t5.*t88)./2.0;
t94 = (t5.*t92)./2.0;
t102 = t69+t75+t76+t82;
t107 = t15.*(t68+t72-t83-t86).*(-1.0./2.0);
t91 = (t3.*t89)./2.0;
t97 = t61+t62+t89;
t98 = t63+t65+t88;
t103 = t68+t72+t78+t81;
t104 = t85+t87+t92;
t99 = (t5.*t98)./2.0;
t100 = (t3.*t97)./2.0;
t106 = (alpha.*t3.*t104)./2.0;
mt1 = [alpha_dot.*(t50+t53+t90-t91+t99-t100+t5.*t104+t3.*(t68+t72-t83-t86)),-beta_dot.*(t50+t53+t99-t100-t3.*t89+t5.*t88+(t5.*t104)./2.0+(t3.*(t68+t72-t83-t86))./2.0)+beta_dot.*(t94+t96+t106+t107),-alpha_dot.*(t50+t53+t90-t91+t99-t100)+alpha_dot.*(t106+t107+t3.*(t68+t72)+(t7.*t89)./2.0+t5.*t92+(t3.*t103)./2.0-(t5.*t102)./2.0+(t15.*t97)./2.0+(alpha.*t3.*t88)./2.0+(alpha.*t3.*t98)./2.0)];
mt2 = [-beta_dot.*(t50+t53+t94+t96+t106+t107+t7.*t89+t15.*t97+alpha.*t3.*t88+alpha.*t3.*t98)-beta_dot.*((t7.*(t68+t72))./2.0+(t15.*t103)./2.0-(alpha.*t3.*t92)./2.0+(alpha.*t3.*t102)./2.0)];
D = reshape([mt1,mt2],2,2);
end
