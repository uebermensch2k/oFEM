function [S,t, area4]=stiffness_matrix_P1_2D(co,el, mat)
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

% read number of coordinates and elements
nC = numel(co(:,1));
%*** First vertex of elements and corresponding edge vectors
X = co(el(:,2),:) - co(el(:,1),:);
Z = co(el(:,3),:) - co(el(:,1),:);
%*** Vector of element areas 4*|T|
area4 = 2*(X(:,1).*Z(:,2)-X(:,2).*Z(:,1));
area4 = area4./mat;

I = el(:,[1 2 3 2 3 3])';
J = el(:,[1 1 1 2 2 3])';
XX = sum(X.*X,2)./area4;
XZ = sum(X.*Z,2)./area4;
ZZ = sum(Z.*Z,2)./area4;
S = [ XX+ZZ-2*XZ, XZ-ZZ, XZ-XX, ZZ, -XZ, XX ]';
S = sparse(I(:),J(:),S(:),nC,nC);
S = S+S'-diag(diag(S));

t=toc;

area4=mat.*area4/4;

end
