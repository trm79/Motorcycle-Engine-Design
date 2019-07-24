%% CONSTANTS

clear all
close all
D = 1776;  %Specified displacement, cc 
N = 4;      %Number of cylinders
r_c = 10;     %Compression ratio
r_s = 1.1;     %S/B stroke to bore ratio 
lambda = .25;    %r/L - .25 to 1/6 typically
R = .287;        %Gas constant, air, kJ/kg*K
k = 1.35;        %Cp/Cv, used in book, average of values at low and high temps
step = .001*pi;    %num iterations = theta*stepsize

%% GEOMETRY
Vd = D * 10^-6/N;     %V displaced, m^3
Vc = Vd/(r_c - 1);    %V clearance, m^3
V1 = Vd + Vc;          %m^3
B = (4*Vd/(pi*r_s))^(1/3);  %Bore diameter, m
S = r_s * B;    %Stroke length, m
r = S/2;        %Crank diameter, m
L = r/lambda;   %Conn rod length

%% INTAKE - constant pressure, temp
T1 = 333;    %starting temperature, book says it's about 30 C above room temp, deg K
P1 = 100;    %starting pressure, kPa
V1 = Vd + Vc;   %m^3
Mm = P1.*V1/(R.*T1);  %Total mass in chamber, kg
theta_suck = [0:step:pi];
P_1(1:length(theta_suck)) = P1;

%% COMPRESSION - isentropic
T2(1) = T1;
P2(1) = P1;
V2(1) = Vc+Vd; %m^3
theta_squeeze = [pi+step:step:2*pi];
r = S/2;
for i = 2:length(theta_squeeze)
V2(i) = Vc + pi*B^2/4*(r*(1-cos(theta_squeeze(i))) + lambda/2*r*sin(theta_squeeze(i))^2);
T2(i) = T2(i-1)*(V2(i-1)/V2(i))^(k-1);       %deg K
P2(i) = P2(i-1)*(V2(i-1)/V2(i))^k;           %kPa
end

T_2 = T2;   %save temp, pressure, volume arrays
P_2 = P2;
V_2 = V2;

P2 = P2(length(theta_squeeze));
T2 = T2(length(theta_squeeze));
V2 = V2(length(theta_squeeze));

%% COMBUSTION - constant volume
T3 = 2552.47; %K - %From CEA
P3 = P2.*(T3./T2); 
V3 = Vc;
 
%% EXPANSION - isentropic
T4(1) = T3;
P4(1) = P3;
V4(1) = Vc;
theta_bang = [2*pi+step:step:3*pi];
for i = 2:length(theta_bang)
V4(i) = Vc + pi*B^2/4*(r*(1-cos(theta_bang(i))) + lambda/2*r*sin(theta_bang(i))^2);
T4(i) = T4(i-1)*(V4(i-1)/V4(i))^(k-1);       %deg K
P4(i) = P4(i-1)*(V4(i-1)/V4(i))^k;           %kPa
end

T_4 = T4;   %save temp, pressure, volume arrays
P_4 = P4;
V_4 = V4;
T4 = T4(length(theta_bang));
P4 = P4(length(theta_bang));
V4 = V4(length(theta_bang));


W34 = Mm* R * (T4 - T3)/(1-k); %kJ
W12 = Mm* R * (T2 - T1)/(1-k); %kJ

Wout = W34+W12;

%% BLOW - exhaust stroke
theta_blow = [3*pi+step:step:4*pi];
P_5(1:length(theta_blow)) = P1;


%% Plot pressure as a function of crank angle
figure
theta = [theta_suck theta_squeeze theta_bang theta_blow];
P_theta = [P_1 P_2 P_4 P_5];
plot(theta, P_theta)
xlabel('Crank Angle (radians)')
ylabel('P (kPa)')
title('Chamber Pressure vs. Crank Angle')


%% Plot T-P Diagram
figure
T_values = [T1,T_2,T3,T_4,T1];
P_values = [P1,P_2,P3,P_4,P1];
plot(T_values, P_values, 'r')
txt1 = '1';
txt2 = '2';
txt3 = '3';
txt4 = '4';
text(T1,P1,txt1)
text(T2,P2,txt2)
text(T3,P3,txt3)
text(T4,P4,txt4)
xlabel('T (K)')
ylabel('P (kPa)')
title('T-P diagram')

%% Plot P-V Diagram
figure
V_values = [V1,V_2,V3,V_4,V1];
P_values = [P1,P_2,P3,P_4,P1];
plot(V_values, P_values, 'r')
txt1 = '1';
txt2 = '2';
txt3 = '3';
txt4 = '4';
text(T1,P1,txt1)
text(T2,P2,txt2)
text(T3,P3,txt3)
text(T4,P4,txt4)
xlabel('V (M^3)')
ylabel('P (kPa)')
title('P-V diagram')

%% Create state table
fprintf('\nState Table\n') 
fprintf('State\tP (kPa)\tV(m^3)\tT(K)\n')
fprintf('1\t%.3f\t%f\t%.3f\n', P1, V1, T1)
fprintf('2\t%.3f\t%f\t%.3f\n', P2, V2, T2)
fprintf('3\t%.3f\t%f\t%.3f *FROM CEA\n', P3, V3, T3)
fprintf('4\t%.3f\t%f\t%.3f\n', P4, V4, T4)

%% Create power output table
fprintf('\n\nWork Output Values\n')
fprintf('RPM\tHP\tMeanPS\tMaxPS\n')
RPM = 5000;
HP = (Wout*RPM/60)/2*N*1.34102;
MeanPS = 2*S*RPM/60;
MaxPS = pi*S*RPM/60;
fprintf('%f\t%f\t%f\t%f\n', RPM, HP, MeanPS, MaxPS)

RPM = 8000;
HP = (Wout*RPM/60)/2*N*1.34102;
MeanPS = 2*S*RPM/60;
MaxPS = pi*S*RPM/60;
fprintf('%f\t%f\t%f\t%f\n\n', RPM, HP, MeanPS, MaxPS)

%% PRESSURE FORCES (X&Y) AS A FUNCTION OF THETA
beta = pi/3;
P_R = P_theta;                          %kPa, Pressure from right cylinder as a function of theta
l = length(P_R);                        %Length of all arrays in this section                   
A = pi*B^2/4;                           %m^2, Bore area
Fpress_max = P3*A;                      %kN, max pressure force, one piston, just cuz it's cool to have

beta_i = 1;                             %Find index of beta in theta array
while theta(beta_i) <= beta
    beta_i = beta_i+1;
end
               
for i = 1:l-beta_i               %kPa, Build P_L array as a function of theta minus beta
    P_L(i) = P_R(beta_i+i);      %Assign everything P_R(1:beta_i) as P_L(l-beta_i:l)
end
for i = (l-beta_i+1):l
    P_L(i) = P_R(l+1-i);
end

for i = 1:l
    F_pressx(i) = A*sin(beta/2)*(P_L(i) - P_R(i));  %kN, Total pressure (L+R) force projected in X dir
    F_pressy(i) = A*cos(beta/2)*(P_L(i) + P_R(i));  %kN, Total Pressure (L+R) force projected in Y dir
    F_pressLx(i) = A*sin(beta/2)*P_L(i);
    F_pressLy(i) = A*cos(beta/2)*P_L(i);
    F_pressRx(i) = A*sin(beta/2)*(- P_R(i));
    F_pressRy(i) = A*cos(beta/2)* P_R(i);
    F_pressL(i) = P_L(i)*A;                         %kN, Pressure from left piston, (dir of piston)
    F_pressR(i) = P_R(i)*A;                         %kN, Pressure from right piston, (dir of piston)
end

figure
plot(theta,F_pressL)
hold on
plot(theta, F_pressR)
plot(theta, F_pressy)
title('Pressure Force')
xlabel('Crank Angle (radians)')
ylabel('Pressure Force (kN)')
legend('Left Piston', 'Right Piston', 'Sum-ydir')


%% RECIPROCATING FORCES (X&Y) AS A FUNCTION OF THETA
Mp = (0.42344+0.18358)/9.81*10^-3;            %kN, reciprocating weight component
Omega = 5000*2*pi/60;       %Radians/sec, angular velocity of crankshaft (N = RPM)    
Z = Mp*r*Omega^2;           
for i = 1:l
    % LEFT CYLINDER
    F_aL1(i) = Z*cos(theta(i)+beta);            %kN, Total Reciprocating force left piston, first order
    F_aL2(i) = Z*lambda*cos(2*(theta(i)+beta)); %kN, Total Reciprocating force left piston, second order
    F_aL1y(i) = cos(beta/2)*F_aL1(i);           %kN, Total Reciprocating force left piston, first order, projected in Y dir
    F_aL1x(i) = -sin(beta/2)*F_aL1(i);           %kN, Total Reciprocating force left piston, first order, projected in Y dir
    F_aL(i) = F_aL1(i)+F_aL2(i);                %kN, Total Reciprocating force first + second order from left piston, (in dir of piston)
    F_aLx(i) = Z*sin(beta)*F_aL(i);
    F_aLy(i) = Z*cos(beta)*F_aL(i);
   
    
    % RIGHT CYLINDER
    F_aR1(i) = Z*cos(theta(i));   
    F_aR2(i) = Z*lambda*cos(2*theta(i));   
    F_aR1y(i) = cos(beta/2)*F_aR1(i);           %kN, Total Reciprocating force left piston, first order, projected in Y dir
    F_aR1x(i) = sin(beta/2)*F_aR1(i);           %kN, Total Reciprocating force left piston, first order, projected in Y dir
    F_aR(i) = F_aR1(i)+F_aR2(i);                %kN, Total Reciprocating force from right piston, (in dir of piston)
    F_aRx(i) = Z*sin(beta)*F_aR(i);
    F_aRy(i) = Z*cos(beta)*F_aR(i);
    
    % SUM X & Y
    F_ax(i) = F_aRx(i) + F_aLx(i);              %kN, Total reciprocating (L+R) force projected in X dir
    F_ay(i) = F_aRy(i) + F_aLy(i);              %kN, Total reciprocating (L+R) force projected in Y dir
end
figure
plot(theta, F_aL1)
hold on
plot(theta, F_aL2)
plot(theta, F_aL)
xlabel('Crank Angle (radians)')
ylabel('Reciprocating Force (kN)')
legend('First Order-left piston', 'Second Order-left piston', 'Sum-left piston')
title('Reciprocating Forces from LEFT Piston - One Throw')

figure
plot(theta, F_aL)
hold on
plot(theta, F_aR)
plot(theta, (F_ay.^2 + F_ax.^2).^(1/2))
xlabel('Crank Angle (radians)')
legend('Sum-left piston', 'Sum-right piston', 'Resultant from both Cylinders')
ylabel('Reciprocating Force (kN)')
title('Resultant Reciprocating Force from both Cylinders - One Throw')

figure
plot(theta, ((F_aL1y+F_aR1y).^2+(F_aL1x+F_aR1x).^2).^0.5)
title('USEFUL')



%% CENTRIPETAL FORCES (X&Y) AS A FUNCTION OF THETA
Mc = (0.815221+0.27216)/9.81*10^-3;                %kN, rotating weight component, for just one throw (1/2crankshaft+Mrot of connecting rod)
Z = Mc*r*Omega^2;
for i = 1:l
    F_cx(i) = Z*sin(beta/2 - theta(i));
    F_cy(i) = Z*cos(beta/2 - theta(i));
    F_c(i) = Z;
end


figure
plot(theta, F_cx)
hold on
plot(theta, F_cy)
plot(theta, F_c)
legend('Centripetal-x', 'Centripetal-y', 'Centripetal-resultant')
xlabel('Crank Angle (radians)')
ylabel('Centripetal Force (kN)')
title('Centripetal Force in Y-dir')

%% COUNTERWEIGHTS

F_cw = 0.4;
for i = 1:l
   %F_nety(i) = F_aR1y(i)+F_aL1y(i)+) - F_cw*cos(theta(i) + beta/2);  %N, counterweight mass * radial location of CG of counterweight
   %F_cw(i) = (F_aR1y(i)+F_aL1y(i)+F_cy(i))*cos(theta(i) + beta/2)/Omega^2;
   %F_cw2(i) = -(F_aR1(i)*cos(theta(i))+F_aL1(i)*cos(theta(i)+beta));
end


figure
% plot(theta, F_cwy)
% hold on
% plot(theta, F_aR1y)
% plot(theta, F_aL1y)
plot(theta,F_nety)
%plot(theta,F_cw2) 
title('AAAHHHHHH')



% %% PLOT X FORCES
% F_x = F_pressx + F_ax + F_cx;
% 
% figure
% plot(theta, F_pressx)
% hold on
% plot(theta, F_ax)
% plot(theta, F_cx)
% plot(theta, F_x)
% legend('Pressure-x', 'Reciprocating-x', 'Centripetal-x', 'Sum-x')
% [F_xmax, Ix] = max(F_x);
% xlabel('Crank Angle (radians)')
% ylabel('Forces in X-dir (kN)')
% title('Total X-dir Forces')
% 
% 
% 
% %% PLOT Y FORCES
% F_y = F_pressy + F_ay + F_cy;
% 
% figure
% plot(theta, F_pressy)
% hold on
% plot(theta, F_ay)
% plot(theta, F_cy)
% plot(theta, F_y)
% [F_ymax, Iy] = max(F_y);
% legend('Pressure-y', 'Reciprocating-y', 'Centripetal-y', 'Sum-y')
% xlabel('Crank Angle (radians)')
% ylabel('Forces in Y-dir (kN)')
% title('Total Y-dir Forces')
% 
% %% FORCES ON CRANKPIN
% for i = 1:l
%     F_crankpin(i) = (F_y(i)^2+F_x(i)^2)^(1/2);
% end
% figure
% plot(theta, F_crankpin)
% 
% [F_crankpinmax, Icrankpinmax] = max(F_crankpin);
% F_crankpinmax
% %theta(Icrankpinmax)
% 
% %% FORCES ON CONNECTING ROD
% sinPhiPlusHalfBeta = R/L*sin(theta)*cos(beta/2)+sqrt(1-(R/L).^2*sin(theta).^2)*sin(beta/2);
% cosPhiPlusHalfBeta = sqrt(1-(R/L).^2*sin(theta).^2)*cos(beta/2)-R/L*sin(theta)*sin(beta/2);
% for i =1 : l
%     Fy_ConRod(i) = (F_aLx(i)+F_aLx(i))*sinPhiPlusHalfBeta(i) + (F_aLy(i)+F_aLy(i))*cosPhiPlusHalfBeta(i) - (F_aRy(i)+F_pressRy(i)-F_cy(i))*cosPhiPlusHalfBeta(i) + (F_cy(i)-F_aRx(i)-F_pressRx(i))*sinPhiPlusHalfBeta(i);
% end
% figure
% plot(theta,Fy_ConRod)
% title('Force on Connecting Rod')
% [Fy_ConRodmax, IConRodmax] = max(Fy_ConRod);
% Fy_ConRodmax
% 
% 
% 
% %% FORCES ON LEFT PISTON PIN
% for i = 1:l
%     F_piston(i) = F_aL(i) + F_pressL(i) + cos(pi/2-beta/2)*(F_cx(i) - sin(beta/2)*(F_aR(i)+F_pressR(i))) + cos(beta/2)*(cos(beta/2)*(F_aR(i)+F_pressR(i))-F_cy(i));
% end
% figure
% plot(theta, F_piston)
% xlabel('Crank Angle (radians)')
% ylabel('Piston Force (kN)')
% title('Sum of all Forces Acting on Left Piston Pin')
% [F_pistonmax, Ipistonmax] = max(F_piston);
% F_pistonmax
% Theta_leftpistonmax = theta(Ipistonmax)         
% 
% 
% 
% 
%     