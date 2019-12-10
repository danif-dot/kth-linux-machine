%%%%%%%%% Milestone 3 %%%%%%%%%%%%

%% small signal parameters

% NMOS 2,5,6

%g_m2= 684E-6;
g_m2= 6.827E-3;
%r_o2= 45.46E3;
r_o2= 4.53E3;
%C_gd2= 2.96E-15;
C_gd2= 29.71E-15;

r_o5= r_o2;

r_o6= r_o2;
C_gd6= C_gd2;

% PMOS 4,7

%r_o4= 91.24E3/3;
r_o4= 9.15E3/3;
%C_gd4= 3.817E-15;
C_gd4= 38.27E-15;

%r_o7= r_o4*3;
r_o7= r_o4*3;
%g_m7= 502.9E-6;
g_m7= 4.98E-3;
C_gd7= C_gd4;
%C_gs7= 8.672E-15;
C_gs7= 86.79E-15;

%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = inv(inv(r_o4)+inv((1+ 1.2*g_m2*r_o5)*r_o2));

R2 = inv(inv(r_o7)+inv(r_o6));

C1= C_gd4 + C_gd2 + C_gs7 + C_gd7*g_m7*R2;

C2 = C_gd6 + C_gd7*(1+ 1/(g_m7*R2));

Adc= g_m2*g_m7*R1*R2;

w_p1= 2*pi*2E6;

Cc= 1/(g_m7*R1*R2*w_p1);

w_p2= g_m7*Cc/(C1*C2+ C2*Cc+ C1*Cc);                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

%Rz=(C2+Cc)/(Cc*g_m7);
%Rz= 1/g_m7 + 1/(w_p2*Cc);
Rz= 1.2*inv(g_m7);

w_z= 1/(Cc*(Rz- 1/g_m7));

%%%%%%%%%%%%%%%%%%%%%%%%
syms w_u
eqn = 1 + (w_u/w_p1)^2 + (w_u/w_p2)^2 + ((w_u)^4)/((w_p1^2)*(w_p2^2)) - (Adc^2)*(1+(w_u/w_z)^2)^2;
sol = solve(eqn, w_u);
W_u = vpa(sol);

th_1 = vpa(atan2(W_u(2), w_p1)*180/pi);

th_2 = vpa(atan2(W_u(2), w_p2)*180/pi);

th_z = vpa(atan2(W_u(2), w_z)*180/pi);

PM = 180 - th_1 - th_2 + th_z