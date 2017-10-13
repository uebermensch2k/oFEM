clear;
close all;

gridPoints=[0.0,0.0;-0.5,-0.5;-0.5,0.5;0.5,-0.5;0.5,0.5];
timePoints=0:0.025:1;
%weights = @(x) 1+0*x(1,:)+0*x(2,:);
weights = @(x) sin(x(1,:))+cos(x(2,:));
coeffs=struct('c_P',0.17,'rho',750,'k',2400,'h',1,'v',[10/3.6;0],'T_a',[1000,100,1000,100],'tMax',0.5);
grid=struct('LL',[-1,-1],'UR',[1,1],'h_max',1/10);
T0=@(x) 300+0*x(1,1,:);
contPos=[-0;-0];
T_b1 = @(x) 0*x(1,1,:);
T_b2 = @(x) -ofem.matrixarray(15000*exp(-10*dot(x,x,1)));
T_b3 = @(x) -ofem.matrixarray(15000*exp(-10*dot(x,x,1)));

% Ambient temperature equals initial temperature
%[T_xt_eq1,mesh]=heatRobin(gridPoints,timePoints,T_b1,T0,coeffs,grid);
%T_xt_eq2=heatRobin(gridPoints,timePoints,T_b2,T0,coeffs,grid);
%T_xt_eq3=heatRobin(gridPoints,timePoints,T_b3,T0,coeffs,grid);

% Ambient temperature greater than initial temperature
%coeffs.T_a=200;
tic
%T_xt_gt1=heatRobin(gridPoints,timePoints,T_b1,T0,coeffs,grid);
%T_xt_gt2=heatRobin(gridPoints,timePoints,T_b2,T0,coeffs,grid);
[T_xt_gt3]=heatRobin(gridPoints,timePoints,T_b3,T0,coeffs,grid);
% [T_xt_gt3,mesh]=heatRobin(gridPoints,timePoints,T_b3,T0,coeffs,mesh);
toc