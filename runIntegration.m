
%------------------INTEGRATION---------------------------------------------

if ~exist('y0', 'var')
    setUpRandomInitialConditions;
end

if exist('simulationTime', 'var')
    tspan = [0 simulationTime];
end

if ~exist('options', 'var')
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
end

%Integrate the System
[t,y] = ode113('nBodyWpar',tspan,y0,options,flag,N,G,Mass);

disp('Integration is done...');