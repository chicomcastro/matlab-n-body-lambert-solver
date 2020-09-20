%% Set up initial conditions

if ~exist("N", "var")
    setUpParameters;
end

% Coordinates given in special system have to be converted to Rectangular.
PI = 3.1415926535898;
c = PI/180;

%% Default values (see setUpSymetricalInitialConditions.m)
Mass=ones(1,N);         % bodies masses

%Initial positions in r, alpha and beta (spherical coordinates)
r0 = ones(1,N);
ar = (0:360/N:360-1)*c;
br = zeros(1,N);

%initial velocities in v, alpha and beta (spherical coordinates)
v0 = 0.6*ones(1,N);
av = 90*ones(1,N)*c;
bv = zeros(1,N)*c;


%% Earth-Sun-Spacechip system
% Setting two firsts bodies to be Sun (red) and Earth (blue)
M_sun = 1.989e30;   % kg
M_earth = 5.972e24; % kg
um = M_sun + M_earth; % kg
Mass(1) = M_sun/um;
Mass(2) = M_earth/um;

UA = 1.496e11;  % m
R_earth_sun = 1*UA;
ud = UA;
r0(1) = 0;
r0(2) = R_earth_sun/ud;
ar(2) = 0;

ut = sqrt(G/G_metric)*sqrt(ud^3/um);  % s

V_earth_sun = sqrt(G_metric*M_sun/R_earth_sun);
v0(1) = 0;
v0(2) = V_earth_sun/ud*ut;

% Spaceship (green)
Mass(3) = 200/um;
r0(3) = r0(2) + 2000e3/ud;
ar(3) = ar(2);
R_mars_sun = 2.2794e11; % m
semieixo = (R_mars_sun + R_earth_sun) / 2;
if ~exist("V_earth_mars", "var")
    V_earth_mars = sqrt(G_metric*M_sun*(1/semieixo+2/R_earth_sun));
end
v0(3) = (V_earth_mars)/ud*ut;
av(3) = av(2);

%% Frame transformation
% Convert positions to cartesean frame
Rec = 0;

for i=1:N
    Rec(3*i-2:3*i,1) = [           ...
        r0(i)*cos(ar(i))*cos(br(i));    ...
        r0(i)*sin(ar(i))*cos(br(i));    ...
        r0(i)*sin(br(i))                ...
    ];
end

A=0;
for i=1:N
    % Convert spherical velocities to rst frame
    A(3*i-2:3*i,1)=[...
        v0(i)*cos(av(i))*cos(bv(i));...
        v0(i)*sin(av(i))*cos(bv(i));...
        v0(i)*sin(bv(i))...
    ];
    % Convert rst velocities to ijk. Each velocity has a different matrix.
    Rec(3*i-2+3*N:3*i+3*N,1) = rstTOijk(Rec(3*i-2:3*i,1))*A(3*i-2:3*i,1);
end

%%
% Atribute initial conditions vector
y0 = Rec;