%N body simulation in inerital frame.

%clear           % clears workspace to get ready for run
%close all

%%
loadData;
setUpParameters;

%%
setUpInitialConditions;

%%
runIntegration;

%%
global shouldPlot
if shouldPlot == 1
    calculateIntegralsOfMotion;
    doPlottings;
end
