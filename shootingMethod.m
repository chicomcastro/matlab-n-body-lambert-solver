clear;
tic
global shouldPlot;
shouldPlot = 0;

%%
loadData;
V_oe_earth = sqrt(G_metric*M_earth/altitude_from_earth);     % [m/s]
V_oe_inertial = V_oe_earth + V_earth_sun;
V_oe_inertial = V_oe_inertial/1e3;  % to [km/s]

t_voo = 2*pi;   % Tempo [rad]
initial_velocity_guess = (V_oe_inertial + 8) * [0 1 0];  % km/s

initial_pos = (R_earth_sun * [1 0 0] + altitude_from_earth * [1 0 0]);  % m
target_pos = R_mars_sun * [-1 0 0]; % m

initial_state.t_voo = t_voo;
initial_state.initial_pos = initial_pos;
initial_state.target_pos = target_pos;
initial_state.v0 = initial_velocity_guess;

x = initial_state;
%%
itr = 1;
while itr <= 1
    disp("---");
    disp("Iteração: " + itr);
    
    itr = itr + 1;
    x1 = x;
    result1 = coust(x1);
    x2 = x;
    x2.v0 = x2.v0 + rand(1,3);
    result2 = coust(x2);
    
    dx1 = result2.x_t2(1) - result1.x_t2(1);
    dx2 = result2.x_t2(2) - result1.x_t2(2);
    dx3 = result2.x_t2(3) - result1.x_t2(3);
    dv1 = result2.v_t1(1) - result1.v_t1(1);
    dv2 = result2.v_t1(2) - result1.v_t1(2);
    dv3 = result2.v_t1(3) - result1.v_t1(3);
    dXdV = [...
        dx1/dv1 dx2/dv1 dx3/dv1; ...
        dx1/dv2 dx2/dv2 dx3/dv2; ...
        dx1/dv3 dx2/dv3 dx3/dv3; ...
    ];
    erro_dist = result1.error;
    delta_velocidade = -1*inv(dXdV)*erro_dist(:);
    x.v0 = x.v0 + delta_velocidade(:)'/1000;
end
toc