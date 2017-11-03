close all;
clear;
clc;

%% USER DEFINED
mesh_file_name = fullfile('coil_plate');
frequency = 50;
Utot      = 1e0;


%% some stuff needed
scl   = 1;%1/50;
omega = 2*pi*frequency;

%% time
% NT    = 1200;
% T_max = 12/frequency;
% T     = linspace(0,T_max,NT);
% dt    = T(2)-T(1);

max_Aphi    =  1e-10;
min_Aphi    = -1e-10;
max_dAphi   =  1e-10;
min_dAphi   = -1e-10;
max_gradphi =  1e-10;
min_gradphi = -1e-10;
max_J       =  1e-10;
min_J       = -1e-10;



%% Mesh import and define boundaries
mesh    = ofem.mesh;
mesh.load_from_inp(mesh_file_name);
mesh.co = mesh.co*scl;
Nc      = size(mesh.co,3);
Ne      = size(mesh.el,1);
% mesh.show;

%% Coils
% coils are indexed by 1 through 4
[~,bd_cp,~,~] = mesh.neumann(1:5);
bd_cp = vertcat(bd_cp{:});

% for plotting coils
Nbd=size(bd_cp,1);
bd_plot=cell(4*Nbd,1);
for j=1:Nbd
	idx=bd_cp(j,:);
    bd_plot((4*j-3):(4*j-2))={squeeze(double(mesh.co(1,1,idx))),squeeze(double(mesh.co(2,1,idx)))};
    bd_plot(4*j)={'r'};
end

coil.max = [ max(mesh.co(:,:,mesh.el(mesh.parts{3,1},:)),[],3); ...
             max(mesh.co(:,:,mesh.el(mesh.parts{3,2},:)),[],3); ...
             max(mesh.co(:,:,mesh.el(mesh.parts{3,3},:)),[],3); ...
             max(mesh.co(:,:,mesh.el(mesh.parts{3,4},:)),[],3) ];

coil.max = reshape(coil.max, 2,[])';

coil.min = [ min(mesh.co(:,:,mesh.el(mesh.parts{3,1},:)),[],3); ...
             min(mesh.co(:,:,mesh.el(mesh.parts{3,2},:)),[],3); ...
             min(mesh.co(:,:,mesh.el(mesh.parts{3,3},:)),[],3); ...
             min(mesh.co(:,:,mesh.el(mesh.parts{3,4},:)),[],3) ];
% copper
coil.min = reshape(coil.min, 2,[])';

coil.kappa  = 5.8e7;
coil.mu     = (1-6.4e-6)*mesh.mu0;

%% Plate
plate.max   = max(mesh.co(:,:,mesh.el(mesh.parts{3,5},:)),[],1);
plate.min   = min(mesh.co(:,:,mesh.el(mesh.parts{3,5},:)),[],1);
% aluminum
plate.kappa = 3.7e7;
plate.mu    = (1+2.2e-5)*mesh.mu0;


%% Space
% air
space.kappa = 0;
space.mu    = mesh.mu0;

%%
kappa      = zeros(Nc,1);
idx        = unique(mesh.el(vertcat(mesh.parts{3,1:4}),:));
kappa(idx) = coil.kappa;
idx        = unique(mesh.el(mesh.parts{3,5},:));
kappa(idx) = plate.kappa;
idx        = unique(mesh.el(mesh.parts{3,6},:));
kappa(idx) = space.kappa;


% mu      = zeros(Nc,1);
% idx     = unique(mesh.el(vertcat(mesh.parts{3,1:4}),:));
% mu(idx) = coil.mu;
% idx     = unique(mesh.el(mesh.parts{3,5},:));
% mu(idx) = plate.mu;
% idx     = unique(mesh.el(mesh.parts{3,6},:));
% mu(idx) = space.mu;


%% finite element
fe = ofem.finiteelement.P1;

qr = ofem.gaussianquadrature(mesh, fe);


%% assemble
opt = struct('S',1,'M',1,'force',@(X) voltage_axi(X,coil,Utot));

eq  = ofem.development.PhiAxi(mesh, fe, qr);
[asm,~,aux] = eq.assemble(opt);


% stiffness
S = (aux.S{1}+aux.S{2}+aux.S{3}+aux.S{4})/coil.mu; % coils
S = S + aux.S{5}/plate.mu; % plate
S = S + aux.S{6}/space.mu; % space

% mass
M = (aux.M{1}+aux.M{2}+aux.M{3}+aux.M{4})*coil.kappa;
M = M+aux.M{5}*plate.kappa;
M = M+aux.M{6}*space.kappa;

% load
b = -(aux.force{1}+aux.force{2}+aux.force{3}+aux.force{4})*coil.kappa;
b = b-aux.force{5}*plate.kappa; % phi is zero here
b = b-aux.force{6}*space.kappa; % phi is zero here


%% system matrix and solution vector
A  = S+1j*omega*M;
u  = sparse(Nc,1);


%% Dirichlet boundary
dirichlet = mesh.dirichlet(6); dirichlet=dirichlet{1};
DOFs      = setdiff(1:Nc,dirichlet);

% u(dirichlet)=0;
% 
% b=b-A*u;


%% solve
u(DOFs) =A(DOFs,DOFs)\b(DOFs);

%% postprocessing
Aphi    =  full(u);
dAphi   =  full(1j*omega*u);
gradphi = -full(b);
% [gradphi_cog,co_cog] = eq.node2cog(gradphi);


%% plot solution over time
NT    = 400;
T_max = 4/frequency;
T     = linspace(0,T_max,NT);

writerObj = VideoWriter('sol_harmonic.avi');
open(writerObj);

figObj=figure;
figObj.Position(3:4)=[4/3, 1]*680;
for i=1:NT
    t=T(i);
    figObj.Name = sprintf('t=%e s, f=%f Hz',t,frequency);

    %% A_{\varphi}
    Aphit = imag(Aphi*exp(1j*omega*t));

 subplot(2,2,1);
    trisurf (mesh.el,double(mesh.co(1,:,:)),double(mesh.co(2,:,:)),Aphit,'FaceColor','interp','EdgeColor','none');
    title('$$A_{\varphi}$$','Interpreter','latex');
    colorbar;
    max_Aphi = max([max_Aphi; Aphit]);
    min_Aphi = min([min_Aphi; Aphit]);
    caxis([min_Aphi, max_Aphi]);
    axis([0 scl 0 scl min_Aphi max_Aphi]);
    hold on
    bd_plot(3:4:4*Nbd)={max(Aphit)*[1;1]};
    plot3(bd_plot{:},'LineWidth',2);
    hold off

    %% \kappa\frac{\partial A_{\varphi}}{\partial t}
%     dAphit = kappa.*imag(1j*omega*Aphi*exp(1j*omega*t));
    dAphit = kappa.*imag(dAphi*exp(1j*omega*t));
    
    subplot(2,2,2);
    trimesh(mesh.el,double(mesh.co(1,:,:)),double(mesh.co(2,:,:)),dAphit,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    title('$$\kappa\frac{\partial A_{\varphi}}{\partial t}$$','Interpreter','latex');
    colorbar;
    max_dAphi = max([max_dAphi; dAphit]);
    min_dAphi = min([min_dAphi; dAphit]);
    caxis([ min_dAphi max_dAphi ]);
    axis([0 scl 0 scl  min_dAphi max_dAphi ]);
%     view (0,90);
    hold on
    bd_plot(3:4:4*Nbd)={max(dAphit)*[1;1]};
    plot3(bd_plot{:},'LineWidth',2);
    hold off
    %% \kappa\nabla\Phi
    gradphit = imag(gradphi*exp(1j*omega*t));

    subplot(2,2,3);
    trisurf (mesh.el,double(mesh.co(1,:,:)),double(mesh.co(2,:,:)),gradphit,'FaceColor','interp','EdgeColor','none');
    title('$$\kappa\nabla\Phi$$','Interpreter','latex');
    colorbar;
    max_gradphi = max([max_gradphi; gradphit]);
    min_gradphi = min([min_gradphi; gradphit]);
    caxis([min_gradphi, max_gradphi]);
%     axis([0 scl 0 scl min_gradphi max_gradphi]);



    %% J
    Jt = -dAphit-gradphit;
    subplot(2,2,4);
    trimesh(mesh.el,double(mesh.co(1,:,:)),double(mesh.co(2,:,:)),Jt,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
    title('$$\vec J = -\kappa\left(\frac{\partial A_{\varphi}}{\partial t}+\nabla\Phi\right)$$','Interpreter','latex');
    view (0,90);
    colorbar;
    max_J = max([max_J; Jt]);
    min_J = min([min_J; Jt]);
    caxis([ min_J max_J ]);

    drawnow;
    writeVideo(writerObj,getframe(figObj));
end

close(writerObj);
