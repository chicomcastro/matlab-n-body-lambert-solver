%--------------------------------------------------------------------------
%----------------------USER SPECIFIED--------------------------------------
%------------------------PARAMETERS----------------------------------------
%--------------------------------------------------------------------------

%As the user you specify these.  The Program does the rest.
%Remember that the system of ODEs is a 6*N system of first order
%scalar ODEs and expect matlab to bog down accordingly.

N=2                   % number of bodies
simulationTime=10     % How long to run the simulation
G=1                   % gravitational constant

%The value of the gravitation constant in the metric system is
%
G_metric=6.67300e-11; % (m3 kg-1 s-2)
%
%which is very small