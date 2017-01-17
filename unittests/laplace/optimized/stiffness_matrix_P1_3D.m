function [S,t]=stiffness_matrix_P1_3D(co,el)
% stiffness_matrix_P1 is the 2D P1 stiffness matrix, i.e.
% it returns the matrix associated with
% $$
% \int_{T_s} \nabla^T \varphi_i(x) \nabla \varphi_j(x)\,dx,
% $$
% with $T_s$ the s-th element and $i,j$ the indices of the nodes.
%
% [S,t]=stiffness_matrix_P1(co,el)
% S is the stiffness matrix, t the time duration for computing S,
% co are the coordinates, el the elements
%

tic

%% read number of coordinates and elements
nC = size(co,1);

%% compute edge vectors
X = co(el(:,2),:) - co(el(:,1),:);
Y = co(el(:,3),:) - co(el(:,1),:);
Z = co(el(:,4),:) - co(el(:,1),:);

%% compute gradients
dphi2 = cross(Y,Z,2);
dphi3 = cross(Z,X,2);
dphi4 = cross(X,Y,2);
dphi1 = -dphi2-dphi3-dphi4;

%% compute areas, 6*|T|=det
area36 = 6*abs(dot(X,dphi2,2));

clear X Y Z;

%% compute stiffness entries
S11 = dot(dphi1,dphi1,2)./area36;
S21 = dot(dphi2,dphi1,2)./area36;
S31 = dot(dphi3,dphi1,2)./area36;
S41 = dot(dphi4,dphi1,2)./area36;

S22 = dot(dphi2,dphi2,2)./area36;
S32 = dot(dphi3,dphi2,2)./area36;
S42 = dot(dphi4,dphi2,2)./area36;

S33 = dot(dphi3,dphi3,2)./area36;
S43 = dot(dphi4,dphi3,2)./area36;

S44 = dot(dphi4,dphi4,2)./area36;

clear dphi1 dphi2 dphi3 dphi4;

%% assembly
I = el(:,[1 2 3 4 2 3 4 3 4 4]);
J = el(:,[1 1 1 1 2 2 2 3 3 4]);
S = [ S11, S21, S31, S41, S22, S32, S42, S33, S43, S44 ];
S = sparse(I(:),J(:),S(:),nC,nC);
S = S+S'-diag(diag(S));


t=toc;

end
