%% Main function
function tests = tests
    tests = functiontests(localfunctions);
end

%% test case setup
function setupOnce(testCase)
    % set tolerance
    testCase.TestData.tol  = 1e-14;
    testCase.TestData.mesh_inp_file = fullfile('meshes','simple_rect.inp');
    
end

%% test loading from inp
function testLoadFromINP(testCase)
    fileName = testCase.TestData.mesh_inp_file;

    fprintf('Loading mesh from inp-file ''%s'' ...',fileName);
    testCase.TestData.mesh = ofem.mesh;
    info = testCase.TestData.mesh.load_from_inp(fileName);
    fprintf(' done. t(read)=%f, t(postproccess)=%f\n',...
        info.time.file_read,...
        info.time.post_proccess);

    % check if valid
end

