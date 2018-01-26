close all
clear

inp_file = 'resonator';

mesh = ofem.mesh;
mesh.load_from_inp(inp_file);

fe = ofem.finiteelement.NE0P;
qr = ofem.gaussianquadrature(mesh,fe);
op = ofem.development.CurlCurl(mesh,fe,qr);

opt.S = 1;
opt.A = 1;
opt.dirichlet{1} = struct('idx',1,'f',0);
opt.dirichlet{2} = struct('idx',2,'f',1);
opt.M = 1;
opt.c = 1;

[asm,info,aux] = op.assemble(opt);

S = asm.S;
M = asm.M;
b = asm.b;
DOFs = asm.DOFs;

u = full(asm.dirichlet);

u(DOFs) = (S(DOFs,DOFs)+M(DOFs,DOFs))\b(DOFs);

uE = u(mesh.el2ed(:,[1:6]));
uE = reshape(uE',6,1,[]);
uE = ofem.matrixarray(uE);
[w,l] = fe.quaddata(3);
phi = fe.phi(l);
[DinvT,detD,D] = mesh.jacobiandata();
cog = mesh.get_cog();

A = D*phi*uE*detD;


















