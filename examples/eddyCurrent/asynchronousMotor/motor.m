% close all;
clear;
clc;

scl=1/100;

mesh_file_name = fullfile('.','motor');

frequency = 50;
omega=2*pi*frequency;
% Utot=1e10;
Utot=100;
T_max = 4*2*pi/omega;

max_Az      =  1e-11;
min_Az      = -1e-11;
max_gradphi =  1e-11;
min_gradphi = -1e-11;
max_J       =  1e-11;
min_J       = -1e-11;
max_dAz     =  1;
min_dAz     = -1;


%% Mesh import and define boundaries
mesh = ofem.mesh;
mesh.load_from_inp(mesh_file_name);
mesh.co = mesh.co*scl;
Nc = size(mesh.co,1);
Ne = size(mesh.el,1);
% mesh.show();

co = double(permute(mesh.co,[3,1,2]));

%% Coils
coil_max_x = [max(co(mesh.el(mesh.parts{3, 7},:),1)); ...
              max(co(mesh.el(mesh.parts{3, 8},:),1)); ...
              max(co(mesh.el(mesh.parts{3, 9},:),1)); ...
              max(co(mesh.el(mesh.parts{3,10},:),1)); ...
              max(co(mesh.el(mesh.parts{3,11},:),1)); ...
              max(co(mesh.el(mesh.parts{3,12},:),1))];
coil_max_y = [max(co(mesh.el(mesh.parts{3, 7},:),2)); ...
              max(co(mesh.el(mesh.parts{3, 8},:),2)); ...
              max(co(mesh.el(mesh.parts{3, 9},:),2)); ...
              max(co(mesh.el(mesh.parts{3,10},:),2)); ...
              max(co(mesh.el(mesh.parts{3,11},:),2)); ...
              max(co(mesh.el(mesh.parts{3,12},:),2))];

coil_min_x = [min(co(mesh.el(mesh.parts{3, 7},:),1)); ...
              min(co(mesh.el(mesh.parts{3, 8},:),1)); ...
              min(co(mesh.el(mesh.parts{3, 9},:),1)); ...
              min(co(mesh.el(mesh.parts{3,10},:),1)); ...
              min(co(mesh.el(mesh.parts{3,11},:),1)); ...
              min(co(mesh.el(mesh.parts{3,12},:),1))];
coil_min_y = [min(co(mesh.el(mesh.parts{3, 7},:),2)); ...
              min(co(mesh.el(mesh.parts{3, 8},:),2)); ...
              min(co(mesh.el(mesh.parts{3, 9},:),2)); ...
              min(co(mesh.el(mesh.parts{3,10},:),2)); ...
              min(co(mesh.el(mesh.parts{3,11},:),2)); ...
              min(co(mesh.el(mesh.parts{3,12},:),2))];
coil.radius = (coil_max_x-coil_min_x)/2;
coil.pos    = [coil_min_x+coil_max_x, coil_min_y+coil_max_y]/2;

% copper
coil.kappa=5.8e7;
coil.mu=(1-6.4e-6)*mesh.mu0;

%% rotor
% aluminum
rotor.kappa=10.02e6;
rotor.mu=1000*mesh.mu0;


%% Space
% air
space.kappa=0;
space.mu=mesh.mu0;


%% finite element
fe = ofem.finiteelement.P1;


%% assemble
op=ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));
opt.S=1;
opt.M=1;
opt.A=[...
       repmat({coil.kappa*coil.mu},1,6), ... % outer coils
       repmat({coil.kappa*coil.mu},1,6), ... % inner coils
       {space.kappa*space.mu}          , ... % air
       {rotor.kappa*rotor.mu}            ... % inner
      ];
opt.dirichlet{1}=struct('idx',1,'f',0);
opt.dirichlet{2}=struct('idx',2,'f',0);
opt.force=@(X) coil.kappa*coil.mu*voltage(X,coil,Utot);
[asm,info,aux]=op.assemble(opt);

% stiffness
S = asm.S;
M = asm.M;

% load
coilforce = [aux.force{7:12}];
% gradphi = -sum(coilforce,2);

kappa = zeros(Nc,1);
idx = unique(mesh.el(vertcat(mesh.parts{3,1:12}),:));
kappa(idx(:)) = coil.kappa;
idx = unique(mesh.el(mesh.parts{3,13},:));
kappa(idx(:)) = space.kappa;
idx = unique(mesh.el(mesh.parts{3,14},:));
kappa(idx(:)) = rotor.kappa;

% u=full(asm.dirichlet);
DOFs = asm.DOFs;


%% solve
phi = linspace(0,2*pi,7); phi = phi(1:6)';

F    = @(t,y) -coilforce(DOFs,:)*sin(omega*t+phi)-S(DOFs,DOFs)*y;

options = odeset('Mass',M(DOFs,DOFs),'OutputFcn',@showprogress_callback);
[T,u]  = ode15s(F,[0,T_max],zeros(numel(DOFs),1),options);

Az = repmat(asm.dirichlet,1,numel(T));
Az(DOFs,:)=u';
clear u;


%% plot solution over time
NT=size(Az,2);

writerObj = VideoWriter('sol_motor.avi');
open(writerObj);

figObj=figure;
figObj.Position(3:4)=[4/3, 1]*680;
opt=[];
for i=2:NT
    t=T(i);
    figObj.Name = sprintf('t=%e s, f=%e Hz',t,frequency);
%     title(sprintf('t=%e s, f=%e Hz',t,frequency));
    %figObj.Title = sprintf('t=%e s, f=%e Hz',t,frequency);

    %% plot vector potential
    Azt=full(Az(:,i));

    subplot(2,3,1);
    trisurf (mesh.el,co(:,1),co(:,2),Azt,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    title('$$A_z$$','Interpreter','latex');
    colorbar;
    max_Az = max([max_Az; Azt]);
    min_Az = min([min_Az; Azt]);
    caxis([min_Az, max_Az]);
    view (0,90);

    %% plot gradient of scalar potential
    gradphit = full(coilforce*sin(omega*t+phi));

    subplot(2,3,2);
    trisurf (mesh.el,co(:,1),co(:,2),gradphit,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    title('$$\kappa \nabla\Phi$$','Interpreter','latex');
    colorbar;
    max_gradphi = max([max_gradphi; gradphit]);
    min_gradphi = min([min_gradphi; gradphit]);
    caxis([min_gradphi, max_gradphi]);

    %% B-field
    Bt = op.gradu(Azt);
    Bt = [Bt(:,2),-Bt(:,1)];

    subplot(2,3,3);
    quiver( co(:,1),co(:,2), Bt(:,1), Bt(:,2), 'LineWidth', 2 );
    title('$$\vec B$$','Interpreter','latex');
    
    %% plot eddy current
    dAzt=kappa.*(full(Azt-Az(:,i-1)));

    subplot(2,3,4);
    trimesh(mesh.el,co(:,1),co(:,2),dAzt,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    title('$$\kappa \frac{\partial A_z}{\partial t}$$','Interpreter','latex');
    colorbar;
    max_dAz = max([max_dAz; dAzt]);
    min_dAz = min([min_dAz; dAzt]);
    caxis([ min_dAz max_dAz ]);

    %% plot current density
    Jt = -(dAzt+gradphit);

    subplot(2,3,5);
    trimesh(mesh.el,co(:,1),co(:,2),Jt,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    title('$$\vec J = -\kappa \left(\frac{\partial A_z}{\partial t}+\nabla\Phi\right)$$','Interpreter','latex');
    colorbar;
    max_J = max([max_J; Jt]);
    min_J = min([min_J; Jt]);
    caxis([ min_J max_J ]);
    
    drawnow;
    writeVideo(writerObj,getframe(figObj));
end

close(writerObj);
