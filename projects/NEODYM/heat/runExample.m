clear;
close all;

gridPoints=[-0.5,-0.5;-0.5,0.5;0.5,-0.5;0.5,0.5];
timePoints=0:0.025:1;
coeffs=struct('c_P',0.17,'rho',750,'k',2400,'h',1,'v',[0/3.6;0],'T_a',100);
grid=struct('LL',[-1,-1],'UR',[1,1],'h_max',1/10);
T0=@(x) 100+0*x(1,1,:);

T_b1 = @(x) 0*x(1,1,:);
T_b2 = @(x) -ofem.matrixarray(15000*exp(-dot(x,x,1)));
T_b3 = @(x) -ofem.matrixarray(15000*exp(-100*dot(x,x,1)));

% Ambient temperature equals initial temperature
[~,T_xt_eq1,mesh]=heatRobin(gridPoints,timePoints,T_b1,T0,coeffs,grid);
%T_xt_eq2=heatRobin(gridPoints,timePoints,T_b2,T0,coeffs,mesh);
%T_xt_eq3=heatRobin(gridPoints,timePoints,T_b3,T0,coeffs,mesh);

% Ambient temperature greater than initial temperature
coeffs.T_a=200;
%T_xt_gt1=heatRobin(gridPoints,timePoints,T_b1,T0,coeffs,mesh);
%T_xt_gt2=heatRobin(gridPoints,timePoints,T_b2,T0,coeffs,mesh);
[~,~,~]=heatRobin(gridPoints,timePoints,T_b3,T0,coeffs,mesh);
