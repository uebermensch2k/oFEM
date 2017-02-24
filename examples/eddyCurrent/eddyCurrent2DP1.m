close all;
clear;
inp_file_name='coax';


%% load mesh
fprintf('Loading mesh ... ');
tic
mesh=ofem.mesh;
mesh.load_from_inp(inp_file_name);
t=toc;
fprintf('done t=%f\n',t);


%% set constants
freq  = 1e6;
mu    = 1*mesh.mu0;
ka_cu = 7e5;
omega = 2*pi*freq;


%% define function space discretization
fe=ofem.finiteelement.P1;
% elem=ofem.finiteelement.P2;
% mesh.extendmesh;


%% define equation type in oFEM
eq=ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));


%% oFEM assmble method
% use the assmble method of the laplace class to assmble all desired
% matrizes and vectors
fprintf('Assembling stiffness matrix (oFEM) ... \n');

% stiffness matrix, mass matrix and right hand side
opt.S = 1;
opt.A = 1/mu;
opt.M = 1;
opt.a = {0,ka_cu,ka_cu};
opt.force = @(x) current(x);
opt.dirichlet = struct('idx',1,'f',0);
[asm,info,~]=eq.assemble(opt);

fprintf('done t=%f\n',info.time2assemble);

%% total assembly and scalar material incorporation
u    = full(asm.dirichlet);
S    = asm.S;
M    = asm.M;
b    = asm.b;
DOFs = asm.DOFs;

%% solve
S       = S-1i*omega*M;
u(DOFs) = S(DOFs,DOFs) \ b(DOFs);


%% plot solution
co = double(permute(mesh.co,[3,1,2]));
figure;
trimesh (mesh.el(:,1:3),co(:,1),co(:,2),real(u),'FaceColor', 'interp', 'EdgeColor', 'black' );
cnt=1;


%% post processing
ka = zeros(size(co,1),1);
ka(unique([mesh.el(mesh.parts{3,2},:); mesh.el(mesh.parts{3,3},:)])) = ka_cu;
    
for t=0:2/(freq*50):2/freq
    
    u_t = u*exp(1i*omega*t);

    eddy_current = 1i*omega*ka.*u;

    stamped_current = zeros( size(co,1), 1 );
    stamped_current(unique(mesh.el(mesh.parts{3,2},:))) = 1e-4/(pi*0.001^2);
    
    total_current = (stamped_current + eddy_current)*exp(1i*omega*t);
    
    B = eq.gradu(real(u_t)); % B = rot A
    B = [B(:,2),-B(:,1)];
    
    B_ABS = sqrt(dot(B,B,2));
    
    
    mesh.export_UCD(fullfile(pwd,'eddy_current'),...
                     strcat('export',num2str(cnt)),...
                     {'Az';1e12*real(u_t);'V'},...
                     {'Magnitude B';B_ABS;'As/m^2'},...
                     {'Total Current';real(total_current);'A'},...
                     {'B';[B, zeros(size(co,1),1)];'As/m^2'});

    cnt = cnt + 1;
end