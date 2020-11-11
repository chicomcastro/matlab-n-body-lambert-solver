%% Set up initial conditions

if ~exist("N", "var")
    setUpParameters;
end

%% Earth-Sun-Spacechip system
% Setting two firsts bodies to be Sun (red) and Earth (blue)
Mass(1) = M_sun/um;
r1 = zeros(1,3);
v1 = zeros(1,3);

% Spaceship (green)
Mass(2) = spaceship_mass/um;
if ~exist("initial_pos", "var")
    initial_pos = (R_earth_sun * [1 0 0] + 300e3 * [1 0 0]);  % m
end
if ~exist("V_exit_earth", "var")
    V_exit_earth = V_hohmann_earth_mars; %km/s
end
r2 = initial_pos/ud;
v2 = (V_exit_earth*1000)/ud*ut;

%%
% Sets initial conditions vector
y0 = [r1, r2, v1, v2];
y0 = y0(:);     % [column vector]