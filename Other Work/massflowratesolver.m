

%Engine Characteristics
rc = 10; %compression ratio
BL = .97; %Bore to stroke ratio
Vd = 870.5; %Volumetric displacement in cubic cm
B = (Vd * 4 * BL / pi) ^ (1/3); %Bore in cm
Vc = Vd / (rc - 1); %Clearance volume (cc) (EQN 2.1)
L = B / BL; %Stroke length (cm)(EQN (2.2)
a = L / 2; %Crank radius (cm)(EQN 2.3)
l = 1.75 * L; %connecting rod length (cm)
theta = 0; %crank angle (initialized value in radians)
s2 = a * cos(theta) + (l^2 - a^2 * sin(theta)^2)^0.5; %distance parameter (EQN2.5)
Vt1 = Vc + (pi / 4) * B^2 * (l + a - s2); %Volume in terms of crank angle (cc) (EQN 2.4)
N = 5000; %Rotational speed (RPM)
Spa = 2 * L * N / 60; %Mean piston speed (cm/s)(Eqn 2.9)

%Chemical/Thermodynamic properties and parameters
FARi = 1/14.7; %Ideal Fuel/Air Ratio
Qlhv = 45.0; %MJ / kg
T0 = 300; %Atmospheric temperature (Kelvin)
P0 = 170275; %Atmospheric pressure (Pa)
pa = 1.9779E-3; %Density of air at T0 and P0 (g/cc)
y1 = 1.403; %ratio of specific heats of intake air
as = 347.44; %Speed of sound through intake air

Mf = 116.775; %molar mass of fuel g/mol
Mo = 31.9988; %molar mass of oxygen g/mol
Mn = 28.013; %molar mass of nitrogen g/mol
   
mA = 12.6; %moles of air per mole of fuel
FAR = Mf / (mA * Mo + mA * 3.76 * Mn); %Actual Fuel/Air Ratio
Xtot = 1 + mA + mA * 3.76; %Total number of moles of reactants
Xf = 1 / Xtot; %Molar mass fraction of fuel
Xo = mA / Xtot; %Molar mass fraction of oxygen
Xn = mA * 3.76 / Xtot; %Molar mass fraction of nitrogen
R = 8.314; %Universal Gas Constant J/molK

R1 = R / (Xf * Mf + Xo * Mo + Xn * Mn); %Gas constant of Air / Fuel mixture kJ/ kg * K
cv1 = 3 * R1; %Coefficient of constant volume (for fuel/air mix)
y = 1.333; %Assume ideal polyatomic gas
cp1 = y * cv1; %Constant pressure specific heat
C1 = P0 * ((Vd + Vc)^y); %Adiabatic pressure/volume constant
Ct1 = P0 * (Vd + Vc) / T0; %Ideal gas equation constant
m = P0 * (Vc + Vd) / (1000 ^ 3 * R1 * T0); %total mass of fuel/air mixture (kg)
mf = m / (1 + (1/FAR)); %Mass of fuel in mixture (kg)
ma = mf / FAR; %Mass of air in mixture (kg)

t = 30 / (pi * N); %Time per half stroke
g = pi * N / 30;
SF = 1.6; %Valve safety factor
Lv = 3 * 0.25 * SF; %Valve lift in cm
Dv = 3; %Valve Diameter in cm
Ac = pi * Lv * Dv / (100 ^ 2); %Valve Curtain Area (m^2)
Cd = 0.509; %Discharge coefficient
mdotchoke = Cd * Ac * P0 / (1000 * ((R1 * T0)^0.5)) * sqrt(y1) * (2 / (y1 + 1)) ^ ((y1 + 1) / (2 * (y1 - 1))); %Choked mass flow rate
Ap = pi * (B / 100) ^ 2 / 4; %Piston crown area m^2
Sp = 2 * (L / 100) * N / 60; %Mean piston speed m/s
Aeavg = 0.6 * Ac * Cd; %Average effective area
Z = Ap * Sp / (Aeavg * as); %Inlet Mach Number
disp(mdotchoke)

for k = 1:180 %Intake stroke (constant pressure and temperature)
   theta(k) = k * pi / 180; %crank angle degrees
   si(k) = a * cos(theta(k)) + (l^2 - a^2 * sin(theta(k))^2)^0.5; %s parameter at theta (EQN2.5)
   Vti(k) = Vc + (pi / 4) * B^2 * (l + a - si(k)); %Volume at theta (EQN2.4)
   mc(k) = P0 / (2 * R1 * T0) * Vti(k) / (1000^3); %Mass in chamber at crank angle theta
   ti(k) = theta(k) * t; %Time at crank angle
   mdot(k) = 0.5 * (pi * P0 * B^2) / (1000^3 * R1 * T0) * (g * (a / 100) * sin(ti(k) * g) + (a / 100)^2 * g * sin(2 * ti(k) * g)); %mass flow rate per time
   vchannel(k) = mdot(k) / (pa * 1000 * (pi * (Dv/100)^2 / 4));
   %vair(k) = mdot(k) / (pa * 1000 * Areaad(k));
end
figure
plot(theta,mdot)
title('Intake Stroke mass flow rate vs Crank Angle')
%Compression stroke - assume adiabatic and reversible
for i = 1:160 %at 180 degrees the piston is at TDC
   thetac(i) = i * pi / 180; %crank angle in radians
   s1(i) = a * cos(thetac(i)) + (l^2 - a^2 * sin(thetac(i))^2)^0.5; %s parameter at theta
   Vt1(i) =  Vd + Vc - (pi / 4) * B^2 * (l + a - s1(i)); %Volume of the fuel mixture at theta (cc)
   Pc(i) = C1 / (Vt1(i)^y); %Pressure at theta (Pascals)
   Tc(i) = Pc(i) * Vt1(i) / Ct1; %Temperature at theta (Kelvin)
end
%Combustion Process
%Constant Volume (pressure and temperature vary)
Qs = Qlhv * mf / m; %Specific Internal Energy Decrease (MJ)
Tad = 2601.55; %Adiabatic Flame Temperature (K)
%Parameters of burned gas mixture
pb = 4.6685; %Density of burned gases (kg/cu m)
R2 = 0.29063; %Gas constant of burned mixture (kJ/kg*K)
cp2 = 1.7329; %Coefficient of constant pressure (for burned gas) kJ / (kg * K)
y2 = 1.2242; %Assume ideal polyatomic gas
 
for ii = 1:20
   theta0(ii) = (ii+160) * pi / 180;
   sc(ii) = a * cos(theta0(ii)) + (l^2 - a^2 * sin(theta0(ii))^2)^0.5; %s parameter at theta
   Vtc(ii) =  Vd + Vc - (pi / 4) * B^2 * (l + a - sc(ii)); %Volume of the fuel mixture at theta (cc)
   Tb(ii) = sind(4.75*ii) * (mf / m) * (Qlhv * 100 / cv1) + Tc(160); %Temperature of gases after combustion
   Pcm(ii) =  (1000^3) * m * R2 * Tb(ii) / Vtc(ii); %Pressure of burned mixture (Pascals)
end

Wc = m * cv1 * (Tb(20) - T0); %Work across entire compression stroke

T1 = Tb(20); %Chamber temperature after combustion (Kelvin)
P2 = m * R2 * (1000 ^ 3) * T1 / Vc; %Chamber pressure after combustion (Pa)
cv2 = cp2 / y2; %Coefficient of constant volume (for burned gas mix)
C2 = P2 * (Vc)^y2; %Adiabatic pressure/volume constant
Ct2 = P2 * (Vc) / T1; %Ideal gas equation constant

%for ij = 1:20
   %theta(ij) = ij * pi / 180;
   %sc2(ij) = a * cos(theta(ij)) + (l^2 - a^2 * sin(theta(ij))^2)^0.5; %s parameter at theta
   %Vtc2(ij) =  Vc + (pi / 4) * B^2 * (l + a - sc2(ij)); %Volume of the fuel mixture at theta (cc)
   %Tb2(ij) = sind(ij) * (mf / m) * (Qlhv * 100 / cv1) + Tb(20); %Temperature of gases after combustion
   %Pcm2(ij) =  (1000^3) * m * R2 * Tb2(ij) / Vtc2(ij); %Pressure of burned mixture (Pascals)
%end

%Expansion stroke - assume adiabatic and reversible
for j = 1:180 %From directly after combustion to 360 which is BDC
   theta1(j) = j * pi / 180; %Crank angle in radians
   s2(j) = a * cos(theta1(j)) + (l^2 - a^2 * sin(theta1(j))^2)^0.5; %s parameter at theta
   Vt2(j) = Vc + (pi / 4) * B^2 * (l + a - s2(j)); %Volume of the fuel mixture at theta (cc)
   Pe(j) = C2 / (Vt2(j)^y2); %Pressure at theta (Pascals)
   Te(j) = Pe(j) * Vt2(j) / Ct2; %Temperature at theta (Kelvin)
   end
We = m * cv2 * (Te(180) - Tb(20)); %Total Expansive Work

for n = 1:180 %Exhaust stroke (constant pressure and temperature)
   theta2(n) = n * pi / 180; %crank angle in radians
   se(n) = a * cos(theta2(n)) + (l^2 - a^2 * sin(theta2(n))^2)^0.5; %s parameter at theta
   Vte(n) =  Vd + Vc - (pi / 4) * B^2 * (l + a - se(n)); %Volume of the gas mixture at theta (cc)
   mce(n) = Pe(160) / (2 * R2 * Te(160)) * Vte(n) / (1000^3); %Mass in chamber at crank angle theta
   te(n) = theta2(n) * t; %Time at crank angle
   mdote(n) = 0.5 * (pi * Pe(160) * B^2) / (1000^3 * R2 * Te(160)) * (g * (a / 100) * sin(te(n) * g) + (a / 100)^2 * g * sin(2 * te(n) * g)); %mass flow rate per time
end
figure
plot(theta2,mdote)
title('Exhuast Stroke mass flow rate vs Crank Angle')

T = [Tc, Tb, Te]; %Temperature across compression and expansion

P = [Pc Pcm Pe]; %Pressure across compression and expansion

V = [Vt1 Vtc Vt2];
thetatotal = [thetac, theta0, theta1];
pmax=max(P)


figure
plot(P,V)
title('Pressure vs Volume')
figure
plot(thetac,Pc)
title('Compression Stroke Process Pressure vs Crank Angle')
figure
plot(theta0,Pcm)
title('Combustion Process Pressure vs Crank Angle')
figure
plot(theta1,Pe)
title('Expansion Stroke Pressure vs Crank Angle')
figure
plot(thetatotal,P)
title('Theta vs Pressure')




