%close all;  clear;
inp_file_name='pipe';

%% load mesh
fprintf('Loading mesh ... ');
tic
mesh=ofem.mesh;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);
% 
% figure
% s=trisurf(mesh.el(:,2:4),mesh.co(:,1),mesh.co(:,2),mesh.co(:,3),zeros(size(mesh.co,1),1));
% xlabel('x');
% ylabel('y');
% zlabel('z');
% set(s,'FaceAlpha',0.1);

%mesh.extendmesh()


%% material data 
% here for steel in GPa!!
% These would be the SI units
% lambda = 1.1538e+11;
% mu     = 7.6923e+10;
% if you had different element blocks, lambda and mu would be vectors of
% approriate size!
lambda = 115.3846;
mu     = 76.9231;


%% define equation type in oFEM

%% define function space discretization
elem=ofem.finiteelement.P1;

qdata= ofem.gaussianquadrature(mesh,elem);
eq=ofem.development.elastic(mesh,elem, qdata);

eq.setmaterial(lambda, mu);


%% oFEM assmble method
% use the assmble method of the elastic class to assmble all desired
% matrizes and vectors. The mass matrix is not implemmted, yet!
fprintf('Assembling stiffness matrix (oFEM) ... \n');

% stiffness matrix and right hand side
opt.S = 1;
%opt.D = 1;
%opt.b = 1;
opt.force{1} = @myload;
opt.dirichlet{1}.data =0;
opt.dirichlet{1}.idx =1;
[asm,info,aux]=eq.assemble(opt);

fprintf('done t=%f\n',info.time2assemble);


%% total assembly and scalar material incorporation

% Material parametrs are scalar values per element part

S = asm.S;


%% get Dirichlet nodes and set free nodes
% Remeber that the ordering is (x,y,z) therfore the node numbers need some
% adjustments to fit the coordinates
 diriNodes = mesh.dirichlet(1);
 dirichlet = diriNodes{1};
 diri_x    = mesh.dim*dirichlet - 2;
 diri_y    = mesh.dim*dirichlet - 1;
 diri_z    = mesh.dim*dirichlet;
 clear diriNodes;
% 
% %% set DOFs (degrees of freedom)
 DOFs = setdiff(1:mesh.dim*size(mesh.co,3),union(union(diri_x, diri_y), diri_z));
% 
% 
% %% solution vector u
 u = zeros(mesh.dim*size(mesh.co,3),1);
% 
% 
% %% Dirichlet boundary data
% u(diri_x) = 0; 
% u(diri_y) = 0;
% u(diri_z) = 0;
% 
% b         = asm.b - S*u;
% 
% 
% %% solve
 u(DOFs) = S(DOFs,DOFs) \ asm.b(DOFs);
%u = S\asm.b;
u=reshape(u,mesh.dim,[])';

%% computing the deformation gradient tensor, stress and strain in each node (F = I + DU/DX)
F       = eq.defGrad(u);
[E,e,s] = eq.StrainStress(u);

t=toc;
fprintf('done t=%f\n',t);
 
%% export 
% mesh.export_UCD(fullfile(pwd,'solution'), 'pipe', {'U',u, 'm'}, ...
%               {'E', E(:,[1,5,9,6,3,2]), '[]'}, {'S', s{1}', 'GPa'});


            
%% plot solution
figure
s=trisurf(mesh.el(:,2:4),double(squeeze(mesh.co(1,:,:))+u(:,1)),double(squeeze(mesh.co(2,:,:))+u(:,2)),...
    double(squeeze(mesh.co(3,:,:))+u(:,3)),zeros(size(mesh.co, 3),1));
set(s,'FaceAlpha',0.5);

