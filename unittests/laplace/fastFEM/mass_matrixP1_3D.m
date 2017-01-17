function M=mass_matrixP1_3D(elements,areas)
NE=size(elements,1); %number of triangles    

%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!  
IP=integration_point_transformation([-0.7236067977, -0.7236067977, -0.7236067977; ...
                                       0.1708203932, -0.7236067977, -0.7236067977; ...
                                      -0.7236067977,  0.1708203932, -0.7236067977; ...
                                      -0.7236067977, -0.7236067977,  0.1708203932]);
weight =[1/4 1/4 1/4 1/4];
phi=shapefun (IP,'P1');
Y_3Dmatrix=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
M_local=zeros(NLB);
for i=1:size(IP,2)
    M_local=M_local+weight(i)*phi(:,i)*phi(:,i)';
end     
M_local_3Dmatrix=astam(areas,reshape(repmat(M_local,1,NE),NLB,NLB,NE));
X_3Dmatrix=permute(Y_3Dmatrix,[2 1 3]);  
M=sparse(X_3Dmatrix(:),Y_3Dmatrix(:),M_local_3Dmatrix(:));