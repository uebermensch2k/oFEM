close all;  clear all;
name='room_with_heater';

%% Variablen

c = 1005.4; % spezific heat capacity [J/(kg*K)] air 1005.4 Steel 470
rho = 1.1839; % Density [kg/m^3] air for 25°C 1.1839 Steel 7.850
lambda = 0.0262 ; %thermal conductivity [W/(m*K)] air 0.0262 Seel 50

%% load mesh
fprintf('Loading mesh ... ');
tic
mesh=ofem.mesh;
mesh.load_from_inp(name);
t=toc;
fprintf('done t=%f\n',t);
co = double(permute(mesh.co,[3,1,2]));
figure
s=trimesh(mesh.el,co(:,1),co(:,2),co(:,3),zeros(size(co,1),1));
set(s,'FaceAlpha',0.1);

%% define function space discretization
fe=ofem.finiteelement.P1;

%% define equation type in oFEM 

eq=ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));

%% oFEM assmble method
% use the assmble method of the laplace class to assmble all desired
% matrizes and vectors
fprintf('Assembling stiffness matrix (oFEM) ... \n');

% only stiffness matrix
% opt.S = 1;
% [asm,info,~]=eq.assemble(opt);

%% stiffness matrix and right hand side
%  opt.S = 1;
%  opt.b = 1;
%  opt.f = @ChargeDensity;
% [asm,info,~]=eq.assemble(opt);

%% stiffness matrix, mass matrix and right hand side
 opt.S = 1;
 opt.M = 1;
 opt.dirichlet{1} = struct('idx',1,'f',290);
 opt.dirichlet{2} = struct('idx',3,'f',300);

[asm,info,~]=eq.assemble(opt);



fprintf('done t=%f\n',info.time2assemble);

%% total assembly and scalar material incorporation

%Matrial parametrs are scalar values per element part
%For convenience the oFEM mesh already has eps0 in SI units stored

S=lambda*asm.S;

M=c*rho*asm.M;
Dofs = asm.DOFs;
u = full(asm.dirichlet);
u_diri = full(asm.dirichlet);


%% Solving

tmin = 0;
tmax = 1000;
dt = 10;

A = S + 1/dt*M;

count=1;
for t=tmin:dt:tmax
    
    % Geschichtsdaten
    u_alt = u;
    
    b = sparse ( size(co,1), 1 );
    
    %    % dirichlet_2 conditions
    % for j = 1 : size(dirichlet_2,1)
    %   g = -1;
    %   b(dirichlet_2(j,:)) = b(dirichlet_2(j,:)) + ...
    %       norm(cross(coneumann(j,3),:)-co(dirichlet_2(j,1),:), ...
    %       co(dirichlet_2(j,2),:)-co(dirichlet_2(j,1),:))) ...
    %       * g/6;
    % end
    
    b = b + 1/dt*M*u_alt - A * u_diri;
    
    u(Dofs) = A(Dofs,Dofs) \ b(Dofs);
   
    t
    mesh.export_UCD(fullfile(pwd,'room_heater'),...
                     strcat('export',num2str(count)),...
                     {'T'; u; 'K'});

    count = count+1;
end


% % Plot the solution
% figure
% s=trisurf(mesh.el(:,2:4),mesh.co(:,1),mesh.co(:,2),mesh.co(:,3),full(u));
% set(s,'FaceAlpha',0.5);
% 
% % Plot the E_Abs
% plot(sqrt(dot(mesh.co,mesh.co,2)),E_abs,'m+');
