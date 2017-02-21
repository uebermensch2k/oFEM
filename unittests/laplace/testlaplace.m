%% Main function
function tests = testlaplace
    tests = functiontests(localfunctions);
end

%% test case setup
function setupOnce(testCase)

    % include optimized
    addpath(fullfile(pwd,'optimized'));

    % include fastFEM
    addpath(fullfile(pwd,'fastFEM'));

    % set tolerance
    testCase.TestData.tol = 1e-12;
    
end

%% Test P1 stiffness and mass 2D
function testP1StiffnessAndMass2D(testCase)
    % load mesh
    mesh = ofem.mesh;
    mesh.hypercube([0,0],[1,1]);
    for i=1:7
        mesh.uniform_refine();
    end
%     mesh.load_from_inp('planarCapacitor');

    % set material
    mat = ones(size(mesh.el, 1),1);

    % load tolerance
    tol  = testCase.TestData.tol;

    co = permute(mesh.co,[3,1,2]);
    el = mesh.el;

    % optimized
    Sopt = stiffness_matrix_P1_2D(co,el,mat);
    Mopt =      mass_matrix_P1_2D(co,el);

    % oFEM
    fe     = ofem.finiteelement.P1;
    solver = ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));
    [asm,~,~]=solver.assemble(struct('S',1,'M',1));

    SoFEM=asm.S;
    MoFEM=asm.M;

    % fastFEM
    [SfastFEM,areas] = stifness_matrixP1_2D(el,co);
    MfastFEM         =     mass_matrixP1_2D(el,areas  );

    verifyEqual(testCase,Sopt,SoFEM,'AbsTol',tol,'Stiffness matrix: Difference between optimized and oFEM exceeds absolute tolerance');
    verifyEqual(testCase,Mopt,MoFEM,'AbsTol',tol,'Mass matrix: Difference between optimized and oFEM exceeds absolute tolerance');

    verifyEqual(testCase,Sopt,SfastFEM,'AbsTol',tol,'Stiffness matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
    verifyEqual(testCase,Mopt,MfastFEM,'AbsTol',tol,'Mass matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
end

% %% Test P2 stiffness and mass 2D
% function testP2StiffnessAndMass2D(testCase)
%     % load mesh
%     mesh = ofem.mesh;
%     mesh.load_from_inp('planarCapacitor');
%     mesh.extendmesh;
% 
%     % set material
%     mat = ones(size(mesh.el, 1),1);
% 
%     % load tolerance
%     tol  = testCase.TestData.tol;
% 
%     % optimized
%     Sopt = stiffness_matrix_P2_2D(mesh.co,mesh.el,mat);
%     Mopt =      mass_matrix_P2_2D(mesh.co,mesh.el);
% 
%     % oFEM
%     elem   = ofem.finiteelement.P2;
%     solver = ofem.laplace(mesh,elem);
%     [asm,~,~]=solver.assemble(struct('S',1,'M',1));
% 
%     SoFEM=asm.S{1};
%     MoFEM=asm.M{1};
%     for i=2:size(asm.S,1)
%         SoFEM=SoFEM+asm.S{i};
%         MoFEM=MoFEM+asm.M{i};
%     end
% 
%     verifyEqual(testCase,Sopt,SoFEM,'AbsTol',tol,'Stiffness matrix: Difference between optimized and oFEM exceeds absolute tolerance');
%     verifyEqual(testCase,Mopt,MoFEM,'AbsTol',tol,'Mass matrix: Difference between optimized and oFEM exceeds absolute tolerance');
% end

%% Test P1 stiffness and mass 3D
function testP1StiffnessAndMass3D(testCase)
    % load mesh
    mesh = ofem.mesh;
    mesh.hypercube([0,0,0],[1,1,1]);
    for i=1:4
        mesh.uniform_refine();
    end
%     mesh = ofem.mesh;
%     mesh.load_from_inp('spaceChargeSphere');

    % set material
%     mat = ones(size(mesh.el, 1),1);

    % load tolerance
    tol  = testCase.TestData.tol;

    co = permute(mesh.co,[3,1,2]);
    el = mesh.el;

    % optimized
    Sopt = stiffness_matrix_P1_3D(co,el);
    Mopt =      mass_matrix_P1_3D(co,el);

    % oFEM
    fe   = ofem.finiteelement.P1;
    solver = ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));
    [asm,~,~]=solver.assemble(struct('S',1,'M',1));

    SoFEM=asm.S;
    MoFEM=asm.M;

    % fastFEM
    [SfastFEM,areas] = stifness_matrixP1_3D(el,co);
    MfastFEM         =     mass_matrixP1_3D(el,areas);

    verifyEqual(testCase,Sopt,SoFEM,'AbsTol',tol,'Stiffness matrix: Difference between optimized and oFEM exceeds absolute tolerance');
    verifyEqual(testCase,Mopt,MoFEM,'AbsTol',tol,'Mass matrix: Difference between optimized and oFEM exceeds absolute tolerance');

    verifyEqual(testCase,Sopt,SfastFEM,'AbsTol',tol,'Stiffness matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
    verifyEqual(testCase,Mopt,MfastFEM,'AbsTol',tol,'Mass matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');

%     sum(sum(abs(SoFEM-SfastFEM)))
%     sum(sum(abs(SoFEM-Sopt)))
%     sum(sum(abs(MoFEM-MfastFEM)))
%     sum(sum(abs(MoFEM-Mopt)))
end

% %% Test P2 mass 3D
% function testP2Mass3D(testCase)
%     % load mesh
%     mesh = ofem.mesh;
%     mesh.load_from_inp('spaceChargeSphere');
%     mesh.extendmesh;
% 
%     % load tolerance
%     tol  = testCase.TestData.tol;
% 
%     % optimized
% %     Sopt = stiffness_matrix_P1_3D(mesh.co,mesh.el);
%     Mopt =      mass_matrix_P2_3D(mesh.co,mesh.el);
% 
%     % oFEM
%     elem   = ofem.finiteelement.P2;
%     solver = ofem.laplace(mesh,elem);
%     [asm,~,~]=solver.assemble(struct('S',0,'M',1));
% 
% %     SoFEM=asm.S{1};
%     MoFEM=asm.M{1};
%     for i=2:size(asm.M,1)
% %         SoFEM=SoFEM+asm.S{i};
%         MoFEM=MoFEM+asm.M{i};
%     end
% 
% %     % fastFEM
% %     [SfastFEM,areas] = stifness_matrixP1_3D(mesh.el,mesh.co);
% %     MfastFEM         =     mass_matrixP1_3D(mesh.el,areas  );
% 
% %     verifyEqual(testCase,Sopt,SoFEM,'AbsTol',tol,'Stiffness matrix: Difference between optimized and oFEM exceeds absolute tolerance');
%     verifyEqual(testCase,Mopt,MoFEM,'AbsTol',tol,'Mass matrix: Difference between optimized and oFEM exceeds absolute tolerance');
% 
% %     verifyEqual(testCase,Sopt,SfastFEM,'AbsTol',tol,'Stiffness matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
% %     verifyEqual(testCase,Mopt,MfastFEM,'AbsTol',tol,'Mass matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
% end
