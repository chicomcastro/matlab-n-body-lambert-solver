%% Set up initial conditions

if ~exist("N", "var")
    setUpParameters;
end

Mass=ones(1,N);  %their masses

%Coordinates given in special system have to be converted to Rectangular.

PI = 3.1415926535898;
c = PI/180;

%Initial positions in r alpha and beta
r0 = ones(1,N);
ar = (0:360/N:360-1)*c;
br = zeros(1,N);

%initial velocities in v alpha and beta.
v0 = 0.6*ones(1,N);
v0(1) = 0;
av = 90*ones(1,N)*c;
bv = zeros(1,N)*c;


%% Earth-Sun system
M_sun = 2e30;   % kg
M_earth = 6e24; % kg
um = M_sun + M_earth; % kg
mu = M_earth/um;
Mass(1) = 1 - mu;
Mass(2) = mu;

UA = 1.496e11;  % m
R_earth_sun = 1*UA;
ud = UA;
r0(1) = 0;
r0(2) = 1;

ut = (ud^3/G_metric/um)^(1/2);  % s

V_earth_sun = sqrt(G_metric*M_sun/R_earth_sun);
v0(1) = 0;
v0(2) = V_earth_sun/ud*ut;

%%
Rec = 0;
%convert positions to cartesean frame

for i=1:N
    Rec(3*i-2:3*i,1) = [           ...
        r0(i)*cos(ar(i))*cos(br(i));    ...
        r0(i)*sin(ar(i))*cos(br(i));    ...
        r0(i)*sin(br(i))                ...
    ];
end

A=0;
for i=1:N
    %convert spherical velocities to rst frame
    A(3*i-2:3*i,1)=[...
        v0(i)*cos(av(i))*cos(bv(i));...
        v0(i)*sin(av(i))*cos(bv(i));...
        v0(i)*sin(bv(i))...
    ];
    %convert rst velocities to ijk.  Each velociety has a different matrix.
    Rec(3*i-2+3*N:3*i+3*N,1) = rstTOijk(Rec(3*i-2:3*i,1))*A(3*i-2:3*i,1);
end

y0 = Rec;