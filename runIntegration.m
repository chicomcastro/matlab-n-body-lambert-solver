
%------------------INTEGRATION---------------------------------------------

if ~exist('y0', 'var')
    setUpRandomInitialConditions;
end

if exist('simulationTime', 'var')
    tspan = [0 simulationTime];
end

if ~exist('options', 'var')
    tstart = tic;
    options=odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@(t,y)myEvent(t,y,tstart));
end

%Integrate the System
[t,y] = ode113(@(t,y)nBodyWpar(t,y,options,flag,N,G,Mass),tspan,y0,options);

function [value, isterminal, direction] = myEvent(t, y, tstart)
    
    value      = toc(tstart) < 1.0;
    isterminal = 1;   % Stop the integration
    direction  = 0;
    
    if value == 0
        disp("Exiting ode113");
    end
end