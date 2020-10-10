%% Set up initial conditions

if ~exist("N", "var")
    setUpParameters;
end

%% Earth-Sun-Spacechip system
% Setting two firsts bodies to be Sun (red) and Earth (blue)
Mass(1) = M_sun/um;
r1 = zeros(1,3);
v1 = zeros(1,3);

Mass(2) = M_earth/um;
r2 = R_earth_sun/ud * [1 0 0];
v2 = V_earth_sun/ud*ut * [0 1 0];

% Spaceship (green)
Mass(3) = spaceship_mass/um;
if ~exist("initial_pos", "var")
    initial_pos = (R_earth_sun * [1 0 0] + 300e3 * [1 0 0]);  % m
end
if ~exist("V_exit_earth", "var")
    semieixo = (R_mars_sun + R_earth_sun) / 2;
    V_exit_earth = sqrt(G_metric*M_sun*(1/semieixo+2/R_earth_sun))/1000*[0 1 0]; %km/s
end
r3 = initial_pos/ud;
v3 = (V_exit_earth*1000)/ud*ut;

%%
% Sets initial conditions vector
y0 = [r1, r2, r3, v1, v2, v3];
y0 = y0(:);     % [column vector]