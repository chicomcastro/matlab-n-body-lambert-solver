%--------------------------------------------------------------------------
%----------------------USER SPECIFIED--------------------------------------
%------------------------PARAMETERS----------------------------------------
%--------------------------------------------------------------------------

%As the user you specify these.  The Program does the rest.
%Remember that the the system of ODEs is a 6*N system of first order
%scalar ODEs and expect matlab to bog down accordingly.


%If G=1 then integration show interesting behavior for 'simulationTime'
%between 2 and 10.  For such values the program run in a resonable amount
%of time for 'N' between 2 and 12.  (Run times often one minute or much
%less, however if the randomly chosen initial conditions pass the system
%near a singularity then the program will run slower even in this range).

%I've run the program with as many as 40 bodies on a home desk top and
%obtained results in an hour.  I've had out of memory errore for more than
%45 bodies (with simulation times over 10).

N=3                  %number of bodies
simulationTime=10     %How long to run the simulation
G=1                  %gravitational constant

%The value of the gravitation constant in the metric system is
%
%G_metric=6.67300 × 10-11 (m3 kg-1 s-2)
%
%which is very small