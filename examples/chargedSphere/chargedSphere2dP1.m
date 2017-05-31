close all;
clear;

inp_file_name='chargedSphere';

%% load mesh
fprintf('Loading mesh ... ');
tic
mesh=ofem.mesh;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);

% figure
% s=mesh.show();
% set(s,'FaceAlpha',0.1);


%% define function space discretization
fe=ofem.finiteelement.P1;


%% define equation type in oFEM 
op=ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));


%% oFEM assmble method
% use the assmble method of the laplace class to assmble all desired
% matrizes and vectors
fprintf('Assembling stiffness matrix (oFEM) ... \n');

% stiffness matrix and right hand side
opt.S         = 1;
opt.D         = 1;
opt.b         = [1;0;0]; 
opt.dirichlet = struct('idx',1,'f',0);
opt.force     = @chargeDensity;
[asm,info,~]  = op.assemble(opt);

fprintf('done t=%f\n',info.time2assemble);


%% total assembly and scalar material incorporation

% Material parametrs are scalar values per element part
% For convenience the oFEM mesh already has eps0 in SI units stored
S    = asm.S;
b    = asm.b;
DOFs = asm.DOFs;


%% solution vector u
% u = zeros(size(mesh.co,3),1);
u = asm.dirichlet;


%% solve
u(DOFs) = S(DOFs,DOFs) \ b(DOFs);


%% compute gradient
E = -op.gradu(u);


%% export
mesh.export_UCD(fullfile(pwd,'chargedSphere'), 'sphere', {'U', u, 'V'}, ...
                {'E', E,'V/m' });

            
%% plot solution
% figure;
% s=tetramesh(mesh.el,double(permute(mesh.co,[3,1,2])),sum(reshape(u(mesh.el),[],mesh.dim+1),2)/(mesh.dim+1));
% set(s,'FaceAlpha',0.5);

% Plot the E_Abs
figure;
plot(squeeze(sqrt(dot(mesh.co,mesh.co,1))),sqrt(dot(E,E,2)),'m+');
