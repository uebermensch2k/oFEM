clear;
close all;

TotalFailed     = 0;
TotalPassed     = 0;
TotalIncomplete = 0;

%% run laplace tests
cd('laplace');
result=runtests('testlaplace.m');
table(result)
TotalPassed     = TotalPassed    +sum(horzcat(result.Passed    ));
TotalFailed     = TotalFailed    +sum(horzcat(result.Failed    ));
TotalIncomplete = TotalIncomplete+sum(horzcat(result.Incomplete));
cd('..');


% %% run elastic tests
% cd('elastic');
% result=runtests('testelastic.m');
% table(result)
% TotalPassed     = TotalPassed    +sum(horzcat(result.Passed    ));
% TotalFailed     = TotalFailed    +sum(horzcat(result.Failed    ));
% TotalIncomplete = TotalIncomplete+sum(horzcat(result.Incomplete));
% cd('..');

% %% run matrixarray tests
% cd('matrixarray');
% result=runtests('testmatrixarray.m');
% table(result)
% TotalPassed     = TotalPassed    +sum(horzcat(result.Passed    ));
% TotalFailed     = TotalFailed    +sum(horzcat(result.Failed    ));
% TotalIncomplete = TotalIncomplete+sum(horzcat(result.Incomplete));
% cd('..');

%% show result
table(TotalPassed,TotalFailed,TotalIncomplete)

