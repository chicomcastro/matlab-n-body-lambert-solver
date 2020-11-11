
%------------------INTEGRATION---------------------------------------------

if ~exist('y0', 'var')
    setUpRandomInitialConditions;
end

if exist('simulationTime', 'var')
    tspan = [0 simulationTime];
end

tstart = tic;
options=odeset('RelTol',1e-8,'AbsTol',1e-8,'Events',@(t,y)myEvent(t,y,tstart));

%Integrate the System
[t,y,te,ye,ie] = ode113(@(t,y)nBodyWpar(t,y,options,flag,N,G,Mass),tspan,y0,options);

function [value, isterminal, direction] = myEvent(t, y, tstart)
    value      = toc(tstart) < 1.0;     % Set integration time limit to 1s
    isterminal = 1;   % Stop the integration
    direction  = 0;
end