%Two body simulation in inerital frame


%clear the workspace
clear;
%output long
format long;


%--------------------------------------------------------------------------
%----------------------NOTES-----------------------------------------------
%--------------------------------------------------------------------------

%Dependencies:
%This program calls 'nBodyWpar.m'

%USE:
%The user adjusts the paremeters below (in the section USER SPECIFIED
%PARAMETERS) to their liking and then saves.  Then type 'twobodyDesign' at
%the matlab command prompt with this program and it's dependencies in the
%search path.

%Computes two body orbits with user specified conic properties.
%The program is set up for working with ellipses, so e is between 0 and 1.
%Program is used in the introduction notes to in the example where the
%elliptic orbits are computed.

%OUTPUT:
%The output is a plot of the two body orbits.  Circle and astric indicate
%the begining and end of the orbit.  These will coincide for elliptic
%orbits if the period calculation is correct.



%PURPOSE:
%Demonstrates that with the methods and formulas derived in the discussion
%of the Kepler problem and the two body problem one can compute conic
%orbits with desired properties.

%--------------------------------------------------------------------------
%----------------------USER SPECIFIED--------------------------------------
%------------------------PARAMETERS----------------------------------------
%--------------------------------------------------------------------------
% DATA

%User Provides
m1=0.75                   %Mass of Primary:      Anything you want
m2=0.25                   %Mass of Secondary:    Anything you want
G=1                       %Gravational Constant: Anything you want
e=0.85                     %eccentricity:         0 < e < 1
d=7.5                     %max distance          Anything you want

%Derived from above or otherwise fixed
M=m1+m2                    %Sum of masses:    don't change
N=2;                       %number of bodies: don't change
Mass=[m1 m2];              %Mass vector:      don't change

%---------------Compute the needed two body initial conditions-------------

%Kepler constants
Ec2=(G*M)^2*(e^2-1)/2
c2=d*G*M*(1-(1+(2*Ec2)/(G*M)^2)^(1/2))
c=sqrt(c2)
thetaDotZero=c/d^2
v_0=d*thetaDotZero

p=c2/(G*M)
a=p/(1-e^2)

%conversion to two body coordinates
x2_0=d/(1+m2/m1)
x1_0=-(m2/m1)*x2_0

x2Dot_0=v_0/(1+m2/m1)
x1Dot_0=-(m2/m1)*x2Dot_0


T=2*pi*a^(3/2)/sqrt(G*M)


%Then the initial conditions are
y0 = [x1_0; 0.0; 0.0;
    x2_0; 0.0; 0.0;
    0.0 ; x1Dot_0; 0.0;
    0.0; x2Dot_0; 0.0];

%integrate
numSteps = 2000;
tf=T;
tspan = linspace(0, tf, numSteps);
%initial positions and velocities
options=odeset('RelTol',1e-13,'AbsTol',1e-22);
[t,y] =  ode113('nBodyWpar',tspan,y0,options,flag,N,G,Mass);

%plotting
hold on

plot(y(numSteps,1),y(numSteps,2),'r*',y(numSteps,4),y(numSteps,5),'b*')
plot(y(1,1),y(1,2),'ro',y(1,4),y(1,5),'bo')
plot(y(:,1),y(:,2), 'r',y(:,4),y(:,5), 'b')




