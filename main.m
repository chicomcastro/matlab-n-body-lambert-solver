%N body simulation in inerital frame.

%clear           % clears workspace to get ready for run
%close all

%%
setUpParameters;

%%
setUpInitialConditions;

%%
runIntegration;

%%
calculateIntegralsOfMotion;

%%
doPlottings;
