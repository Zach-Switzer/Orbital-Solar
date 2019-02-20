% GW-Sat Available Power Calculations

% Created: 2/19 - Zach Switzer

% Credits & References:
%   - "An Evaluation of CubeSat Orbital Decay" by D. Oltrogge, K. Leveque
%   (https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1144&context=smallsat)
%   - [1]"How To Do a 3D Circle in MATLAB" by Kye Taylor 
%   (https://www.mathworks.com/matlabcentral/answers/79296-how-to-do-a-3d-circle-in-matlab)
%   - "Trace Marker Along a Line" by MathWorks 
%   (https://www.mathworks.com/help/matlab/creating_plots/trace-marker-along-line.html)
%   - [2]"Earth-sized Sphere With Topography" by Will Campbell 
%   (https://www.mathworks.com/matlabcentral/fileexchange/27123-earth-sized-sphere-with-topography)

% Assumptions:
%   - CubeSat is tumbling in space such that it is difficult to 
%     determine the exact attitude of the CubeSat and it is assumed that 
%     each side / combination of sides has an equal probability of being
%     exposed to the sun at any given time during orbit (while outside of
%     the dark zone of orbit.
%   - Only up to 3 sides can be visible at a time and only adjacent sides
%     can be visible at the same time

clear all; 
close all;
clc;

%% Constants - All Units SI

cell_area=2.66e-3; % [m^2]
s1=6*cell_area; % [m^2]
s2=6*cell_area; % [m^2]
s3=6*cell_area; % [m^2]
s4=5*cell_area; % [m^2]
s5=2*cell_area; % [m^2]
s6=0*cell_area; % [m^2]

% Find all the side combos:
CSA=[s1 s2 s3 s4 s5 s6 ...
    (0.5*(s1+s2)) (0.5*(s1+s4)) (0.5*(s2+s3)) (0.5*(s3+s4))...
    (0.5*(s1+s2+s5)) (0.5*(s1+s4+s5)) (0.5*(s2+s3+s5)) (0.5*(s3+s4+s5))...
    (0.5*(s1+s2+s6)) (0.5*(s1+s4+s6)) (0.5*(s2+s3+s6)) (0.5*(s3+s4+s6))]; 

CSA_avg=mean(CSA);
fprintf('The average exposed area during tumbling is:\n%f [m^2] \n\n', CSA_avg)

%% Orbit Info

Alt=[100 200 300 350 400 500 600]; % Altitude [km]
Mu=3.986e5; % Standard Gravitational Parameter of Earth [km^3/s^2]
R=6371; % Radius of Earth [km]

V=zeros(1,length(Alt)); % Creating an Array of Zeros for the For Loop
T=zeros(1,length(Alt));
for i=1:length(Alt)
    
    V(:,i)=(sqrt(Mu./(Alt(i)+6371)))*1000; % Velocity Based on Orbit Altitude Assuming Circular Orbit [m/s]
    T(:,i)=2*pi*sqrt(((R+Alt(i)).^3)/Mu)/3600; % Orbital period (circular) [h]

end

% Finding the orbit time in light and dark
rev=24./T; % Number of revolutions / day
theta=acosd(R./(Alt+R)); % Angle between the radius of the earth and the orbit [deg]
phi=(180-2*theta)*pi/180; % Angle of shadow behind the Earth [rad]
arc_len=(6371+Alt).*phi; % Length of "Darkness" arc [km]
c=2*pi.*(6371+Alt); % Length of total orbial [km]
per_dark=(arc_len./c).*100; % Percentage of orbit in Darkness [%]
per_light=((c-arc_len)./c).*100; % Percentage of orbit in Light [%]
time_dark=(per_dark./100).*T; % Time in darkeness per orbit [h]
time_light=(per_light./100).*T; % Time in light per orbit [h]

% Assume: Guaranteed 14-16 revolutions per day but don't know what will
% happen with the fractional value because we don't know the orbital
% starting position every day. 
%   - The only way to get an accurate number for this would be to compute
%     the new position after each day, and take into account orbital decay
%     due to drag. 
%   ^ The method above will be added on at a later time. 

rev_guar=floor(rev); % Guaranteed number of orbits per day
light_day=time_light.*rev_guar; % Time in light / day [h]
dark_day=time_dark.*rev_guar; % Time in dark / day [h]
unknown_d=24-(light_day+dark_day); % Unknown time in light &/or dark / day [h]

%% Power Generation

% Constants
BOL=29.5; % Beginning of life SolAero cell efficiency [%]
EOL_good=26.55; % End of life SolAero cell efficiency best-case [%]
EOL_bad=25.67; % End of life SolAero cell efficiency worst-case [%]
SI=[1353 1361 1367]; % Solar Irradiance Values found online [W/m^2]
sun_power=mean(SI); % Mean Solar Irradiance [W/m^2]

% Calculations
power_BOL=sun_power*CSA_avg*(BOL/100); % Power @ BOL [W]
power_EOL_good=sun_power*CSA_avg*(EOL_good/100); % Power @ EOL best-case [W]
power_EOL_bad=sun_power*CSA_avg*(EOL_bad/100); % Power @ EOL worst-case [W]
energy_BOL=power_BOL*time_light; % Energy gained / orbit BOL [Wh]
energy_EOL_good=power_EOL_good*time_light; % Energy gained / orbit EOL best-case [Wh]
energy_EOL_bad=power_EOL_bad*time_light; % Energy gained / orbit EOL worst-case [Wh]
power_tot=[power_BOL, power_EOL_good, power_EOL_bad]; 
power_avg=mean(power_tot); % Average power gained / orbit [W]
energy_tot=[energy_BOL, energy_EOL_good, energy_EOL_bad];
energy_avg=mean(energy_tot); % Average energy gained / orbit [Wh]

fprintf('The average power gained in orbit is:\n%f [W] \n\n', power_avg)
fprintf('The average energy gained in orbit is:\n%f [Wh] \n\n', energy_avg)

% Daily #'s
power_avg_day=power_avg.*rev_guar; % Average power gained / day [W]
energy_avg_day=power_avg.*light_day; % Average energy gained / day [Wh]
%fprintf('The average energy gained per day is:\n%f [Wh] \n\n', energy_avg_day)

%% Plotting Data

% 1-D Plots:
% Time in Light & Darkness / Orbit 
figure;
plot(Alt,time_dark, '-r*') 
title('Time in Light & Dark Per Orbit'); 
xlabel('Altitude [km]');
ylabel('Time [h]'); 
legend('show','Location','northwest');
hold on;
plot(Alt,time_light,'-b*') 
legend('Darkness', 'Light') 
%axis([100 600 65 90]);

% Time in Light & Darkness / Day
figure;
plot(Alt,dark_day, '-r*') 
title('Time in Light & Dark Per Day'); 
xlabel('Altitude [km]');
ylabel('Time [h]'); 
legend('show','Location','northwest');
hold on;
plot(Alt,light_day,'-b*') 
legend('Darkness', 'Light') 

% Daily Energy & Power
figure;
title('Energy Gained Per Day'); 
xlabel('Altitude [km]');
ylabel('Watt-Hours [Wh]'); 
legend('show','Location','northwest');
hold on;
plot(Alt,energy_avg_day,'-r*') 
legend('Average Energy') 
%axis([100 600 65 90]);

%% Orbital / Animation 

% Uses earth_sphere.m file 
% => Make sure this is in the same folder as this file when running

% Create the orbit => level circle 
theta=linspace(0,2*pi,150); %values to run through
rad=6371+600; %radius
xo=rad*cos(theta); %x-values of circle eqn
yo=rad*sin(theta); %y-values of circle eqn
zo=0*theta;

% [1] Create the orbit => rotated circle
% Original points, original plane
x = rad*cos(theta);
y = rad*sin(theta);
z = 0*theta;
pnts = [x;y;z];
% unit normal for original plane
n0 = [0;0;1]; 
n0 = n0/norm(n0);

% unit normal for plane to rotate into 
% plane is orthogonal to n1... given by equation
% n1(1)*x + n1(2)*y + n1(3)*z = 0
%n1 = [1;0;100]; 
n1 = [2;-5;4];
n1 = n1/norm(n1); 
% theta is the angle between normals
c = dot(n0,n1) / ( norm(n0)*norm(n1) ); % cos(theta)
s = sqrt(1-c*c);                        % sin(theta)
u = cross(n0,n1) / ( norm(n0)*norm(n1) ); % rotation axis...
u = u/norm(u); % ... as unit vector
C = 1-c;
% the rotation matrix
rot = [u(1)^2*C+c, u(1)*u(2)*C-u(3)*s, u(1)*u(3)*C+u(2)*s
    u(2)*u(1)*C+u(3)*s, u(2)^2*C+c, u(2)*u(3)*C-u(1)*s
    u(3)*u(1)*C-u(2)*s, u(3)*u(2)*C+u(1)*s, u(3)^2*C+c];
% Rotated points
newPnts = rot*pnts;

n2 = [-3;5;4];
n2 = n2/norm(n2); 
% theta is the angle between normals
c2 = dot(n0,n2) / ( norm(n0)*norm(n2) ); % cos(theta)
s2 = sqrt(1-c2*c2);                        % sin(theta)
u2 = cross(n0,n2) / ( norm(n0)*norm(n2) ); % rotation axis...
u2 = u2/norm(u2); % ... as unit vector
C2 = 1-c2;
% the rotation matrix
rot2 = [u2(1)^2*C2+c2, u2(1)*u2(2)*C2-u2(3)*s2, u2(1)*u2(3)*C2+u2(2)*s2
    u2(2)*u2(1)*C2+u2(3)*s2, u2(2)^2*C2+c2, u2(2)*u2(3)*C2-u2(1)*s2
    u2(3)*u2(1)*C2-u2(2)*s2, u2(3)*u2(2)*C2+u2(1)*s2, u2(3)^2*C2+c2];
% Rotated points
newPnts2 = rot2*pnts;

% [2] 3-D Plots & Animations
figure;
earth_sphere('km');
hold on;
plot3(newPnts(1,:),newPnts(2,:),newPnts(3,:),'k', 'LineWidth', 1.0)
plot3(xo,yo,zo,'-r', 'LineWidth', 1.0)
title('Orbit Visualization')
plot3(newPnts2(1,:),newPnts2(2,:),newPnts2(3,:),'m', 'LineWidth', 1.0)
p = plot(xo(1),yo(1),'sr','MarkerFaceColor','r','MarkerSize',10);

for k = 2:length(theta)
    p.XData = xo(k);
    p.YData = yo(k);
    drawnow
end

% follow this format or it'll look weird when you plot it
p2 = plot(newPnts(1,:),newPnts(2,:),'sk','MarkerFaceColor','k','MarkerSize',10);
for f = 2:length(theta)
    p2.XData = newPnts(1,f);
    p2.YData = newPnts(2,f);
    p2.ZData = newPnts(3,f);
    drawnow
end

p3 = plot(newPnts2(1,:),newPnts2(2,:),'sm','MarkerFaceColor','m','MarkerSize',10);
for g = 2:length(theta)
    p3.XData = newPnts2(1,g);
    p3.YData = newPnts2(2,g);
    p3.ZData = newPnts2(3,g);
    drawnow
end