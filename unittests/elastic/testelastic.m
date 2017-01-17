%% Main function
function tests = testelastic
    tests = functiontests(localfunctions);
end

%% test case setup
function setupOnce(testCase)

    % include fastFEM
    addpath(fullfile(pwd,'fastFEM'));

    % set tolerance
    testCase.TestData.tol = 1e-12;

    % lame constants
    testCase.TestData.lambda = 115.3846;
    testCase.TestData.mu     = 76.9231;
    
end

%% Test P1 stiffness and mass 2D
function testP1StiffnessAndMass2D(testCase)
    % load mesh
    mesh = ofem.mesh;
    mesh.load_from_inp('myplate2');

    % set lame constants
    lambda = testCase.TestData.lambda;
    mu     = testCase.TestData.mu    ;

    % load tolerance
    tol  = testCase.TestData.tol;

    % oFEM
    elem   = ofem.finiteelement.P1;
    solver = ofem.elastic(mesh,elem);
    solver.setmaterial(lambda, mu);
    [asm,~,~]=solver.assemble(struct('S',1,'M',1));

    SoFEM=asm.S{1};
    MoFEM=asm.M{1};
    for i=2:size(asm.S,1)
        SoFEM=SoFEM+asm.S{i};
        MoFEM=MoFEM+asm.M{i};
    end

    % fastFEM
    [SfastFEM,areas] = stifness_matrixP1_2D_elasticity(mesh.el,mesh.co,lambda,mu);
    MfastFEM         =     mass_matrixP1_2D_elasticity(mesh.el,areas);

    verifyEqual(testCase,SoFEM,SfastFEM,'AbsTol',tol,'Stiffness matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
    verifyEqual(testCase,MoFEM,MfastFEM,'AbsTol',tol,'Mass matrix: Difference between fastFEM and oFEM exceeds absolute tolerance');
end

