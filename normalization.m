% function [G, um, ud, ut] = normalizeUnits(Mass, Distance)
% Makes unit normalization for n body problem

% Chooses gravitational constant as 1
G = 1;  % gravitational constant normalized

% Set unit of mass as the sum of all masses involved
% In this case, it's specificied for Three Body Restricted Problem (3BRP)
um = M_sun + M_venus;	% [kg]

% Set unit of distance as the distance between two mainly bodies
% In this case, it's specificied for Sun-Earth system
ud = UA;	% [m]

% Set unit of time as the expected to have G = 1
ut = sqrt(G/G_metric)*sqrt(ud^3/um);	% [s]
