%% Main function
function tests = testmatrixarray
    tests = functiontests(localfunctions);
end

%% test case setup
function setupOnce(testCase)
    % set tolerance
    testCase.TestData.tol = 1e-14;
    
end

% %% Test 2D operations
% function testOperations2D(testCase)
%     % tolerance
%     tol = testCase.TestData.tol;
% 
%     % create data
%     A=rand(2,1,10);
%     B=ofem.matrixarray(A);
% 
%     % invoke tests
%     verifyEqual(testCase,raw(rot(B)),[-A(2,1,:); A(1,1,:)],'AbsTol',tol,'matrixarray.rot failed');
% end
% 
% %% Test 3D operations
% function testOperations3D(testCase)
% %     % tolerance
% %     tol = testCase.TestData.tol;
% % 
% %     % create data
% %     A=rand(3,1,10);
% %     B=ofem.matrixarray(A);
% end
% 
% %% Test general operations
% function testOperations(testCase)
%     % tolerance
%     tol = testCase.TestData.tol;
% 
%     N = 10;
%     M = 20;
%     K=10;
% 
%     % create data
%     % matrices
%     A=ofem.matrixarray(rand(N,M,K));
%     B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(M,N));
%     D=ofem.matrixarray(rand(N,M));
% 
%     % vectors
%     u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
%     w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));
% 
%     % scalars
%     a=ofem.matrixarray(rand(1,1));
%     b=ofem.matrixarray(rand(1,1,K));
% 
%     %% test matrix + matrix
%     
% 
% 
%     A=rand(N,M,10);
%     B=ofem.matrixarray(A);
% 
%     C=rand(M,N);
%     D=rand(N,M);
%     E=rand(M,N,10);
%     a=rand(1,1);
%     v=rand(M,1);
%     w=rand(1,N);
% 
%     %% test vector x matrix
% %     clear res2;
%     res1=raw(B');
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)';
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
% 
%     %% test matrix x matrix
%     clear res2;
%     res1=raw(B*C);
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)*C;
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'mtimes(matrixarray,matrix) failed');
% 
%     %% test matrix x matrix
%     clear res2;
%     res1=raw(B*E);
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)*E(:,:,i);
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'mtimes(matrixarray,matrixarray) failed');
% 
%     %% test matrix x matrix
%     clear res2;
%     res1=raw(E*B);
%     for i=1:10
%         res2(:,:,i)=E(:,:,i)*A(:,:,i);
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'mtimes(matrixarray,matrixarray) failed');
% 
%     %% test matrix x vector
%     clear res2;
%     res1=raw(B*v);
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)*v;
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'mtimes(matrixarray,vector) failed');
% 
%     %% test matrix x scalar
%     clear res2;
%     res1=raw(B*a);
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)*a;
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'mtimes(matrixarray,vector) failed');
% 
%     %% test vector x matrix
%     clear res2;
%     res1=raw(w*B);
%     for i=1:10
%         res2(:,:,i)=w*A(:,:,i);
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'mtimes(vector,matrixarray) failed');
% 
%     %% test matrix + matrix
%     clear res2;
%     res1=raw(B+D);
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)+D;
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'plus(matrixarray,matrix) failed');
% 
%     %% test matrix - matrix
%     clear res2;
%     res1=raw(B-D);
%     for i=1:10
%         res2(:,:,i)=A(:,:,i)-D;
%     end
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'plus(matrixarray,matrix) failed');
% 
%     %% index operations
%     res1=B(:,:,end);
%     res2=A(:,:,end);
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'B(:,:,1) failed');
%     res1=raw(B(end,:,:));
%     res2=A(end,:,:);
%     verifyEqual(testCase,res1,res2,'AbsTol',tol,'B(end,:,:) failed');
% end

%% test addition
function testOperationsPlus(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    N = 10;
    M = 20;
    K = 10;

    %% create data
    % matrices
    A=ofem.matrixarray(rand(N,M,K));
    B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(N,M));
    D=ofem.matrixarray(rand(M,N));

    % vectors
    u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
    w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));

    % scalars
    a=ofem.matrixarray(rand(1,1,K));
    b=ofem.matrixarray(rand(1,1));

    %% matrix + matrix
    clear res2;
    res1=double(A+B');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))+double(B(:,:,i))';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A+D');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))+double(D)';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix + scalar
    clear res2;
    res1=double(A+a);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))+a(i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A+b);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))+b;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + vector
    clear res2;
    res1=double(u+u);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))+double(u(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u+w);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))+double(w);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + scalar
    clear res2;
    res1=double(u+a);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))+double(a(i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u+b);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))+double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% scalar + scalar
    clear res2;
    res1=double(a+a);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))+double(a(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(a+b);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))+double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end

%% test subtraction
function testOperationsMinus(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    N = 10;
    M = 20;
    K = 10;

    %% create data
    % matrices
    A=ofem.matrixarray(rand(N,M,K));
    B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(N,M));
    D=ofem.matrixarray(rand(M,N));

    % vectors
    u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
    w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));

    % scalars
    a=ofem.matrixarray(rand(1,1,K));
    b=ofem.matrixarray(rand(1,1));

    %% matrix + matrix
    clear res2;
    res1=double(A-B');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))-double(B(:,:,i))';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A-D');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))-double(D)';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix + scalar
    clear res2;
    res1=double(A-a);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))-a(i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A-b);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))-b;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + vector
    clear res2;
    res1=double(u-u);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))-double(u(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u-w);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))-double(w);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + scalar
    clear res2;
    res1=double(u+a);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))+double(a(i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u-b);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))-double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% scalar + scalar
    clear res2;
    res1=double(a-a);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))-double(a(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(a-b);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))-double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end


%% test .*
function testOperationsTimes(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    N = 10;
    M = 20;
    K = 10;

    %% create data
    % matrices
    A=ofem.matrixarray(rand(N,M,K));
    B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(N,M));
    D=ofem.matrixarray(rand(M,N));

    % vectors
    u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
    w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));

    % scalars
    a=ofem.matrixarray(rand(1,1,K));
    b=ofem.matrixarray(rand(1,1));

    %% matrix + matrix
    clear res2;
    res1=double(A.*B');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).*double(B(:,:,i))';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A.*D');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).*double(D)';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix + scalar
    clear res2;
    res1=double(A.*a);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).*a(i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A.*b);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).*b;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + vector
    clear res2;
    res1=double(u.*u);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).*double(u(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u.*w);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).*double(w);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + scalar
    clear res2;
    res1=double(u.*a);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).*double(a(i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u.*b);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).*double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% scalar + scalar
    clear res2;
    res1=double(a.*a);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i)).*double(a(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(a.*b);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i)).*double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end


%% test *
function testOperationsmTimes(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    N = 10;
    M = 20;
    K = 10;

    %% create data
    % matrices
    A=ofem.matrixarray(rand(N,M,K));
    B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(N,M));
    D=ofem.matrixarray(rand(M,N));

    % vectors
    u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
    w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));

    % scalars
    a=ofem.matrixarray(rand(1,1,K));
    b=ofem.matrixarray(rand(1,1));

    %% matrix * matrix
    clear res2;
    res1=double(A*B);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))*double(B(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A*D);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))*double(D);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix * vector
    clear res2;
    res1=double(A*u);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))*u(:,1,i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A*w);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))*w;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix * scalar
    clear res2;
    res1=double(A*a);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))*a(i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A*b);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))*b;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector * vector
    clear res2;
    res1=double(u*u');
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))*double(u(:,:,i)');
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u*w');
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))*double(w');
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector * scalar
    clear res2;
    res1=double(u*a);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))*double(a(i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u*b);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))*double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% scalar * scalar
    clear res2;
    res1=double(a*a);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))*double(a(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(a*b);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))*double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end


%% test ./
function testOperationsrDivide(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    N = 10;
    M = 20;
    K = 10;

    %% create data
    % matrices
    A=ofem.matrixarray(rand(N,M,K));
    B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(N,M));
    D=ofem.matrixarray(rand(M,N));

    % vectors
    u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
    w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));

    % scalars
    a=ofem.matrixarray(rand(1,1,K));
    b=ofem.matrixarray(rand(1,1));

    %% matrix + matrix
    clear res2;
    res1=double(A./B');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))./double(B(:,:,i))';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A./D');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))./double(D)';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix + scalar
    clear res2;
    res1=double(A./a);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))./a(i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A./b);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i))./b;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + vector
    clear res2;
    res1=double(u./u);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))./double(u(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u./w);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))./double(w);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + scalar
    clear res2;
    res1=double(u./a);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))./double(a(i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u./b);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i))./double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% scalar + scalar
    clear res2;
    res1=double(a./a);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))./double(a(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(a./b);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i))./double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end


%% test .\
function testOperationslDivide(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    N = 10;
    M = 20;
    K = 10;

    %% create data
    % matrices
    A=ofem.matrixarray(rand(N,M,K));
    B=ofem.matrixarray(rand(M,N,K));
%     C=ofem.matrixarray(rand(N,M));
    D=ofem.matrixarray(rand(M,N));

    % vectors
    u=ofem.matrixarray(rand(M,1,K));
%     v=ofem.matrixarray(rand(1,N,K));
    w=ofem.matrixarray(rand(M,1));
%     x=ofem.matrixarray(rand(1,N));

    % scalars
    a=ofem.matrixarray(rand(1,1,K));
    b=ofem.matrixarray(rand(1,1));

    %% matrix + matrix
    clear res2;
    res1=double(A.\B');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).\double(B(:,:,i))';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A.\D');
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).\double(D)';
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% matrix + scalar
    clear res2;
    res1=double(A.\a);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).\a(i);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(A.\b);
    for i=1:K
        res2(:,:,i)=double(A(:,:,i)).\b;
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + vector
    clear res2;
    res1=double(u.\u);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).\double(u(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u.\w);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).\double(w);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% vector + scalar
    clear res2;
    res1=double(u.\a);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).\double(a(i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(u.\b);
    for i=1:K
        res2(:,:,i)=double(u(:,:,i)).\double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% scalar + scalar
    clear res2;
    res1=double(a.\a);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i)).\double(a(:,:,i));
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    clear res2;
    res1=double(a.\b);
    for i=1:K
        res2(:,:,i)=double(a(:,:,i)).\double(b);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end


%% test subtraction
function testOperationsRot(testCase)
    % tolerance
    tol = testCase.TestData.tol;

    K = 10;

    %% create data
    % vectors
    u=ofem.matrixarray(rand(2,1,K));
    v=ofem.matrixarray(rand(1,2,K));

    %% rot(u)
    clear res2;
    res1=double(rot(u));
    res2=zeros(2,1,K);
    for i=1:K
        res2(:,:,i)=double([-u(2,:,i);u(1,:,i)]);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');

    %% rot(v)
    clear res2;
    res1=double(rot(v));
    for i=1:K
        res2(:,:,i)=double([-v(:,2,i),v(:,1,i)]);
    end
    verifyEqual(testCase,res1,res2,'AbsTol',tol,'transpose failed');
end
