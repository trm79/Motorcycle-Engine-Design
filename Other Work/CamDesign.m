function CamDesign_Intake
clear all

N = 5000;       %RPM
B = 0.103;     %bore diameter, m
theta = [0:0.001:4*pi];
phi(1) = 0;
maxlift = 0.016;

for i = 2:length(theta)
    phi(i) = phi(i-1) + 2*pi/length(theta); %phi = angle
end

start = pi/2
finish = pi
midpoint = (start + finish) / 2; % midpoint change from rise to fall
dist = finish - start;

%Lift
for i = 1:length(theta)
    lift(i) = (maxlift/2)*(((1-cos((phi(i)*(1/(dist/(2*pi))))))));
    if phi(i) < start
        lift(i) = 0;
    end
    if phi(i) > finish
      lift(i) = 0;
    end
end
    
plot(phi, lift)
xlabel('Angle (radians)')
ylabel('Lift of Valve (m)')

%Cycloidal
beta = midpoint;

for i = 1:length(theta)
    s(i) = (maxlift/2)*((1-cos(pi*phi(i)/beta))-0.25*(1-cos(2*pi*phi(i)/beta)));
    v(i) = ((pi*maxlift)/(beta*2))*(sin(pi*phi(i)/beta)-0.5*sin(2*pi*phi(i)/beta));
    a(i) = (((pi^2)*maxlift)/((beta^2)*2))*(cos(pi*phi(i)/beta)-cos(2*pi*phi(i)/beta));
    j(i) = ((-(pi^3)*maxlift)/((beta^3)*2))*(sin(pi*phi(i)/beta)-2*sin(2*pi*phi(i)/beta));
    %if phi(i) < start
    %    s(i) = 0;
     %   v(i) = 0;
      %  a(i) = 0;
       % j(i) = 0;
    %end
    %if phi(i) > finish
    %    s(i) = 0;
     %   v(i) = 0;
      %  a(i) = 0;
       % j(i) = 0;
    %end
end

figure
plot(phi, s)
hold on
plot(phi, v)
plot(phi, a)
plot(phi, j)
legend('s (m)', 'v (m/rad)', 'a (m/rad^2)','j(m/rad^3)')
title('S, V, A, J - Cycloidal')
xlabel('Angle (rad)')
axis([0 2*pi -0.1 0.1])

figure
plot(phi, lift)
hold on
plot(phi, s)
legend('Lift', 'Cycloidal')
ylabel('Lift (m)')
xlabel('Angle (rad)')

%Cam Profile
r_base = 0.05; %base radius, m

%Flat Follower
for i = 1 : length(theta)
    x(i) = (r_base+lift(i))*cos(phi(i));
    y(i) = (r_base+lift(i))*sin(phi(i));
end

%Base Circle Plot
for i = 1 : length(phi)
    xbase(i) = r_base*cos(phi(i));
    ybase(i) = r_base*sin(phi(i));
end

figure
plot(x,y)
hold on
plot(xbase,ybase)
%axis([-0.05 0.05 -0.05 0.05])
axis square
title('Cam Profile')




