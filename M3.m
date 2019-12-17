%%%%%%%%% Milestone 3 %%%%%%%%%%%%
clear all
close all
%% small signal parameters

%=========NMOS 1,2,5,6=========

%===1,2===%
%W/L=100u/150n  ID=400u
g_m2= 6.827E-3;
r_o2= 4.53E3;
C_gd2= 29.71E-15;

%W/L=10u/150n  ID=40u
%g_m2= 684E-6;
%r_o2= 45.46E3;
%C_gd2= 2.96E-15;

%W/L=10u/500n  ID=200u
%g_m2= 1.467E-3;
%r_o2= 23.96E3;
%C_gd2= 3.21E-15;

%===6===%
%W/L=100u/150n  ID=400u
r_o6= 4.53E3;
C_gd6= 29.71E-15;

%W/L=10u/500n  ID=200u
%r_o6= 23.96E3;
%C_gd6= 3.21E-15;

%===5===%
%W/L=200u/150n  ID=800u
r_o5= 2.259E3;

%W/L=10u/500n  ID=400u
%r_o5= 9.797E3;

%=========PMOS 3,4,7===========
%===3,4===%
%W/L=100u/150n  ID=400u
r_o4= 9.15E3;
C_gd4= 38.27E-15;
g_m4 = 4.98E-3; 

%W/L=  ID=
%r_o4= 91.24E3/3;
%C_gd4= ;

%===7===%
%W/L=100u/150n  ID=400u
r_o7= 9.15E3;
g_m7= 4.98E-3;
C_gd7= 38.27E-15;
C_gs7= 86.79E-15;

%W/L=10u/150n  ID=40u
%r_o7= 91.24E3;
%g_m7= 502.9E-6;
%C_gd7= 3.817E-15;
%C_gs7= 8.672E-15;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%R1 = inv(inv(r_o4)+inv((1+ 1.2*g_m2*r_o5)*r_o2));
%R1 = inv(inv(r_o2)+inv(r_o4));
r_d= inv(inv(1/g_m4)+inv(r_o4));
R1= (2*r_o2+r_d)*r_o4*inv((1+g_m4*r_d)*r_o4+2*r_o2+r_d);

R2 = inv(inv(r_o7)+inv(r_o6));
%R2 = 6E3;
CL = 1E-12;

C1= C_gd4 + C_gd2 + C_gs7 + C_gd7*g_m7*R2;

C2 = C_gd6 + C_gd7*(1+ 1/(g_m7*R2));

A1= g_m2*R1;
A2= g_m7*R2;
Adc= A1*A2;
Adc_dB= 20*log10(Adc);

%% previous approach %%
%w_p1= 2*pi*2E6;
%f_p1 = w_p1/(2*pi);

%Cc= 1/(g_m7*R1*R2*w_p1);

%w_p2= g_m7*Cc/(C1*C2+ C2*Cc+ C1*Cc);                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
%f_p2 = w_p2/(2*pi);

%Rz=(C2+Cc)/(Cc*g_m7);
%Rz= 1/g_m7 + 1/(w_p2*Cc);  %Rz equal to gm_7
%Rz= 2*inv(g_m7);            %Rz greater than g_m7 => zero on LHP

%w_z= 1/(Cc*(Rz - 1/g_m7));
%f_z = w_z/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%

%% new approach %%
Cc = 1.3E-12;

Rz =  (CL+Cc)*inv(g_m7*Cc);

w_p1 = inv(R1*((1+g_m7*R2)*(Cc+C_gd7)+C1)+R2*(Cc+C_gd7+CL));
f_p1 = w_p1*inv(2*pi);

w_p2 = inv(w_p1)*inv(R1*R2*((Cc+C_gd7)*C1+(Cc+C_gd7)*CL+C1*CL));
f_p2 = w_p2*inv(2*pi);

w_z = 1/(Cc*(Rz - 1/g_m7));
f_z = w_z*inv(2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%
% syms w_u
% %eqn = 1 + (w_u/w_p1)^2 + (w_u/w_p2)^2 + ((w_u)^4)/((w_p1^2)*(w_p2^2)) - (Adc^2)*(1+(w_u/w_z)^2)^2;
% eqn = 1 + (w_u/w_p1)^2 + (w_u/w_p2)^2 + ((w_u)^4)/((w_p1^2)*(w_p2^2)) - (Adc^2)*(1+(w_u/w_z)^2)^2;
% sol = solve(eqn, w_u);
% W_u = double(max(real(vpa(sol))));
% f_u = W_u/(2*pi);

% th_1 = atan2(W_u, w_p1)*180/pi;
% 
% th_2 = atan2(W_u, w_p2)*180/pi;
% 
% th_z = atan2(W_u, w_z)*180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%

%% calculating w_u %%
%pol = [1, 0, w_p1^2 + w_p2^2, 0, ((w_p1*w_p2)^2)*(1- Adc^2)];
pol = [w_z^2, 0, w_z^2 *(w_p1^2 + w_p2^2)- (Adc*w_p1*w_p2)^2, 0, (w_z*w_p1*w_p2)^2*(1- Adc^2)];
X = roots(pol);
w_u = max(real(X));
f_u = w_u*inv(2*pi);

th_1 = atan2(w_u, w_p1)*180/pi;

th_2 = atan2(w_u, w_p2)*180/pi;

th_z = atan2(w_u, w_z)*180/pi;

PM = 180 - th_1 - th_2 + th_z;
%%%%%%%%%%%%%%%%%%%%%%

GBW = Adc*w_u;


