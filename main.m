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
    %% Quantities needed for the constants of motion loops
    A=size(y);                 %vector containing (rows,columns)
    timeMax=A(1,1);             %number of rows of y is number of time steps
    M=0;                       %initialize M
    for i=1:N
        M=M+Mass(1,i);   %sum of the masses is needed to
                         %compute the constantsclo of motion
    end
    %%
    %calculateIntegralsOfMotion;
    %%
    doPlottings;
end
