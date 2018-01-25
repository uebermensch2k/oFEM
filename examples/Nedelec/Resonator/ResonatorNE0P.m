inp_file = 'resonator';

mesh = ofem.mesh;
mesh.load_from_inp(inp_file);

fe = ofem.finiteelement.NE0P;
qr = ofem.gaussianquadrature(mesh,fe);
op = ofem.development.CurlCurl(mesh,fe,qr);

opt.S = 1;
opt.A = 1;
opt.dirichlet{1} = struct('idx',1,'f',0);
opt.M = 1;
opt.c = mesh.mu0;

[asm,info,aux] = op.assemble(opt);
