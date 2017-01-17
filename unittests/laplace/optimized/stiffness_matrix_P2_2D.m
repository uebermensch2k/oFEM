function [S,t]=stiffness_matrix_P2_2D(co,el, mat)
% stiffness_matrix_P2 is the 2D P2 stiffness matrix, i.e.
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

% read number of coordinates and elements
nC = numel(co(:,1));
%*** First vertex of elements and corresponding edge vectors
X = co(el(:,2),:) - co(el(:,1),:);
Y = co(el(:,3),:) - co(el(:,2),:);
Z = co(el(:,1),:) - co(el(:,3),:);
%*** Vector of element areas 2*|T|=det(.)
area2 = (-X(:,1).*Z(:,2)+X(:,2).*Z(:,1));
area2=area2.*mat;

I = el(:,[1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6])';
J = el(:,[1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6])';
XX = sum(X.*X,2)./area2;
XY = sum(X.*Y,2)./area2;
XZ = sum(X.*Z,2)./area2;
YY = sum(Y.*Y,2)./area2;
YZ = sum(Y.*Z,2)./area2;
ZZ = sum(Z.*Z,2)./area2;
OO = zeros(size(XX));
S = [ YY/2,  -YZ/6,  -XY/6,  2*YZ    /3,    OO      , 2* XY    /3, ...
     -YZ/6,   ZZ/2,  -XZ/6,  2*YZ    /3, 2* XZ    /3,    OO      , ...
     -XY/6,  -XZ/6,   XX/2,    OO      , 2* XZ    /3, 2* XY    /3, ...
    2*YZ/3, 2*YZ/3,   OO  , 4*(XX-YZ)/3, 4* XY    /3, 4* XZ    /3, ...
      OO  , 2*XZ/3, 2*XZ/3, 4* XY    /3, 4*(YY-XZ)/3, 4* YZ    /3, ...
    2*XY/3,   OO  , 2*XY/3, 4* XZ    /3, 4* YZ    /3, 4*(ZZ-XY)/3]';

S = sparse(I(:),J(:),S(:),nC,nC);
% S = S+S'-diag(diag(S));

t=toc;

end
