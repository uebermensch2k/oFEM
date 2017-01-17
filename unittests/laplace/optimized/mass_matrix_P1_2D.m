function [M,t]=mass_matrix_P1_2D(co,el)
% mass_matrix_P1_2D is the 2D P1 mass matrix, i.e.
% it returns the matrix associated with
% $$
% \int_{T_s} \varphi_i(x) \varphi_j(x)\,dx,
% $$
% with $T_s$ the s-th element and $i,j$ the indices of the nodes.
%
% [M,t]=mass_matrix_P1_2D(co,el)
% M is the mass matrix, t the time duration for computing M,
% co are the coordinates, el the elements
%

tic

% read number of coordinates and elements
nC = numel(co(:,1));
%*** First vertex of elements and corresponding edge vectors
X = co(el(:,2),:) - co(el(:,1),:);
Z = co(el(:,1),:) - co(el(:,3),:);
%*** Vector of element areas 2*|T|==det(.)
a = (-X(:,1).*Z(:,2)+X(:,2).*Z(:,1))'/24;

I = el(:,[1 2 3 2 3 3])';
J = el(:,[1 1 1 2 2 3])';
M = [2;1;1;2;1;2]*a;
M = sparse(I(:),J(:),M(:),nC,nC);
M = M+M'-diag(diag(M));

t=toc;

end
