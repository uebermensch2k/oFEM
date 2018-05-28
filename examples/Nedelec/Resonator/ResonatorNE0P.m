close all
clear

inp_file = 'resonator';

freq = 1e9;

mesh = ofem.mesh;
mesh.load_from_inp(inp_file);

fe = ofem.finiteelement.NE0P;
qr = ofem.gaussianquadrature(mesh,fe);
op = ofem.development.CurlCurl(mesh,fe,qr);

opt.S = 1;
opt.A = 1/mesh.mu0;
opt.dirichlet{1} = struct('idx',1,'f',0);
opt.dirichlet{2} = struct('idx',2,'f',0);
opt.M = 1;
opt.c = mesh.eps0*freq;

[asm,info,aux] = op.assemble(opt);

S = asm.S;
M = asm.M;
b = asm.b;
DOFs = asm.DOFs;

u = full(asm.dirichlet);
u(DOFs) = rand(size(DOFs,1),1);

u(DOFs) = (S(DOFs,DOFs)+M(DOFs,DOFs))\b(DOFs);

uE = u(mesh.el2ed(:,[1:6]));
uE = reshape(uE',6,1,[]);
uE = ofem.matrixarray(uE);
[w,l] = fe.quaddata(3);
phi = fe.phi(l);
[DinvT,detD,D] = mesh.jacobiandata();
cog = mesh.get_cog();

A = D*phi*uE*(1/detD);
Aabs = sqrt(A(1,1,:).^2+A(2,1,:).^2+A(3,1,:).^2);

x = cog(1,1,:);
x = double(x(:));
y = cog(2,1,:);
y = double(y(:));
z = cog(3,1,:);
z = double(z(:));
Aabs = double(Aabs(:));

scatter3(x,y,z,1,(Aabs))
















