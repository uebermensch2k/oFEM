close all;
clear all;
inp_file_name='coax';


%% set constants
freq  = 1e6;
mu_r  = 1;
ka_cu = 7e5;
omega = 2*pi*freq;


%% load mesh
fprintf('Loading mesh ... ');
tic
mesh=ofem.mesh;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);


%% define function space discretization
elem=ofem.finiteelement.P1;
% elem=ofem.finiteelement.P2;
% mesh.extendmesh;


%% define equation type in oFEM
eq=ofem.laplace(mesh,elem);


%% oFEM assmble method
% use the assmble method of the laplace class to assmble all desired
% matrizes and vectors
fprintf('Assembling stiffness matrix (oFEM) ... \n');

% stiffness matrix, mass matrix and right hand side
opt.S = 1;
opt.M = 1;
opt.b = 1;
opt.load = @current;
[asm,info,~]=eq.assemble(opt);

fprintf('done t=%f\n',info.time2assemble);

%% total assembly and scalar material incorporation

% Material parametrs are scalar values per element part
% For convenience the oFEM mesh already has eps0 in SI units stored

mu = mesh.mu0*mu_r;

S=(asm.S{1}+asm.S{2}+asm.S{3}+asm.S{4})/mesh.mu0;
M=ka_cu*asm.M{2}+ka_cu*asm.M{3};

%% get Dirichlet nodes and set free nodes
diriNodes = mesh.dirichlet(1);
dirichlet = diriNodes{1};
clear diriNodes;


%% set DOFs (degrees of freedom)
k = setdiff(1:size(mesh.co,1),dirichlet); %degrees of freedom


%% solution vector u
u = zeros(size(mesh.co,1),1);


%% Dirichlet boundary data
u(dirichlet) = 0;
% b            = asm.b - S*u;
b            = asm.b;


%% solve
S       = S-1i*omega*M;
u(k) = S(k,k) \ b(k);


%% plot solution
figure;
trimesh ( mesh.el(:,[1 2 3]), mesh.co(:,1), mesh.co(:,2), real(u), 'FaceColor', 'interp', 'EdgeColor', 'black' );
cnt=1;


%% post processing
for t=0:2/(freq*50):2/freq
    
    u_t = u*exp(1i*omega*t);
    
    ka = zeros( size(mesh.co,1), 1 );
    ka(unique([mesh.el(mesh.parts{3,2},:); mesh.el(mesh.parts{3,3},:)])) = ka_cu;
    
    eddy_current = 1i*omega*ka.*u;

    stamped_current = zeros( size(mesh.co,1), 1 );
    stamped_current(unique(mesh.el(mesh.parts{3,2},:))) = 1e-4/(pi*0.001^2);
    
    total_current = (stamped_current + eddy_current)*exp(1i*omega*t);
    
    B = eq.gradu(real(u)); % B = rot A
    B = [B(:,2),-B(:,1)];
    
    B_ABS = sqrt(dot(B,B,2));
    
    
    mesh.export_UCD(fullfile(pwd,'eddy_current'),...
                     strcat('export',num2str(cnt)),...
                     {'Az';1e12*real(u_t);'V'},...
                     {'Magnitude B';B_ABS;'As/m^2'},...
                     {'Total Current';real(total_current);'A'},...
                     {'B';[B, zeros(size(mesh.co,1),1)];'As/m^2'});

    cnt = cnt + 1;
end