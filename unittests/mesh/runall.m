clear;
close all;

TotalFailed     = 0;
TotalPassed     = 0;
TotalIncomplete = 0;

%% run matrixarray tests
result=runtests('tests.m');
table(result)
TotalPassed     = TotalPassed    +sum(horzcat(result.Passed    ));
TotalFailed     = TotalFailed    +sum(horzcat(result.Failed    ));
TotalIncomplete = TotalIncomplete+sum(horzcat(result.Incomplete));

%% show result
table(TotalPassed,TotalFailed,TotalIncomplete)

