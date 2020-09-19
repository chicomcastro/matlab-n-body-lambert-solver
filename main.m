%N body simulation in inerital frame.

clear           % clears workspace to get ready for run

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
