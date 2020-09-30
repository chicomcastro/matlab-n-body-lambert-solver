
M_sun = 1.989e30;   % kg
M_earth = 5.972e24; % kg
spaceship_mass = 200; % kg

UA = 1.496e11;  % m
R_earth_sun = 1*UA;
R_mars_sun = 2.2794e11; % m

altitude_from_earth = 2000e3;   % m


%The value of the gravitation constant in the metric system is
%
G_metric=6.67300e-11; % (m3 kg-1 s-2)
%
%which is very small

V_earth_sun = sqrt(G_metric*M_sun/R_earth_sun); % m/s

% Coordinates given in special system have to be converted to Rectangular.
PI = 3.1415926535898;
rad = PI/180;

% Rotation matrix
Rx = @(t) [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = @(t) [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz = @(t) [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];

normalization;