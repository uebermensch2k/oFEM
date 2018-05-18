close all;
clear;
inp_file_name='planar_trelis';
msh_file_name='planar_gmsh';


%% set constants
U_1    = -1;
U_2    =  1;

eps1_r = 1; % Left
eps2_r = 2; % Right
eps3_r = 1; % Room


%% load mesh
fprintf('Loading trelis file ... \n');
tic
meshtrelis=ofem.mesh;
meshtrelis.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);

fprintf('Loading gmsh file ... \n');
tic
meshgmsh=ofem.mesh;
meshgmsh.load_from_msh(msh_file_name);
t=toc;
fprintf('done t=%f\n',t);


fe=ofem.finiteelement.P1;
qrtrelis = ofem.gaussianquadrature(meshtrelis,fe);
optrelis=ofem.elliptic(meshtrelis,fe,qrtrelis);

qrgmsh = ofem.gaussianquadrature(meshgmsh,fe);
opgmsh = ofem.elliptic(meshgmsh,fe,qrgmsh);

opt.S            = 1;
opt.A            = {eps1_r*meshtrelis.eps0,eps2_r*meshtrelis.eps0,eps3_r*meshtrelis.eps0};
opt.dirichlet{1} = struct('idx',1,'f',U_1);
%opt.dirichlet{2} = struct('idx',2,'f',U_2);
opt.neumann{1}   = struct('idx',2,'f',1e-8);
%opt.neumann{2}   = struct('idx',1,'f',U_1);
[asmtrelis,~,auxtrelis]     = optrelis.assemble(opt);
[asmgmsh,~,auxgmsh]     = opgmsh.assemble(opt);


S    = asmtrelis.S;
DOFs = asmtrelis.DOFs;
b    = asmtrelis.b;

utrelis = full(asmtrelis.dirichlet);

tic
utrelis(DOFs) = S(DOFs,DOFs) \ b(DOFs);
t=toc;

S    = asmgmsh.S;
DOFs = asmgmsh.DOFs;
b    = asmgmsh.b;

ugmsh = full(asmgmsh.dirichlet);

tic
ugmsh(DOFs) = S(DOFs,DOFs) \ b(DOFs);
t=toc;



%% plot solution
figure(1);
h=optrelis.plot(utrelis);
set(h,'FaceColor', 'interp', 'EdgeColor', 'black');

figure(2);
h=opgmsh.plot(ugmsh);
set(h,'FaceColor', 'interp', 'EdgeColor', 'black');

% figure(3)
% spy(asmtrelis.S);
% 
% figure(4)
% spy(asmgmsh.S);

