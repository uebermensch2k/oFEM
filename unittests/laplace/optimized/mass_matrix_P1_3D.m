function [M,t]=mass_matrix_P1_3D(co,el)
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
% edges
e12 = co(el(:,2),:) - co(el(:,1),:);
e13 = co(el(:,3),:) - co(el(:,1),:);
e14 = co(el(:,4),:) - co(el(:,1),:);
% det
detD = dot(e12,cross(e13,e14,2),2)';


I = el(:,[1 2 3 4 2 3 4 3 4 4])';
J = el(:,[1 1 1 1 2 2 2 3 3 4])';
M = [2/120;1/120;1/120;1/120;2/120;1/120;1/120;2/120;1/120;2/120]*detD;
M = sparse(I(:),J(:),M(:),nC,nC);
M = M+M'-diag(diag(M));

t=toc;

end
