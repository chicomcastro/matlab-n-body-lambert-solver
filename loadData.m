% script loadData
% Loads all physics quantities for running simulation, as planet masses
% and positions, spaceship requisites, gravitation constant, and utils as
% each 3D rotation matrices, SOI calculator, pi, and rad2degree const

%The value of the gravitation constant in the metric system is
%
G_metric=6.67300e-11; % (m3 kg-1 s-2)
%
%which is very small

% Masses
M_sun = 1.989e30;   % kg
M_earth = 5.972e24; % kg
M_mars = 6.39e23;   % kg
spaceship_mass = 200; % kg

% Distances
UA = 1.496e11;  % m
R_earth_sun = 1*UA;
R_mars_sun = 2.2794e11; % m

% Altitude of parking orbit before exit Earth
altitude_from_earth = 400e3;   % m

% Earth orbit velocity
V_earth_sun = sqrt(G_metric*M_sun/R_earth_sun);             % [m/s]

% Parking orbit from Earth
V_oe_earth = sqrt(G_metric*M_earth/altitude_from_earth);    % [m/s]
V_oe_inertial = V_oe_earth + V_earth_sun;                   % [m/s]
V_oe_inertial = V_oe_inertial/1e3;                          % [km/s]

% Hohmann transfer velocity from Earth to Mars
semieixo = (R_mars_sun + R_earth_sun) / 2;
V_exit_earth = sqrt(G_metric*M_sun*(1/semieixo+2/R_earth_sun))/1000*[0 1 0]; %km/s

% Coordinates given in special system have to be converted to Rectangular.
PI = 3.1415926535898;
rad = PI/180;

% Rotation matrix
Rx = @(t) [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = @(t) [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz = @(t) [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

% Sphere of influence
SOI = @(a,m,M) a*(m/M)^(2/5);
SOI_mars = SOI(R_mars_sun, M_mars, M_sun);
SOI_earth = SOI(R_earth_sun, M_earth, M_sun);

% Makes unit normalization for N body problem
normalization;