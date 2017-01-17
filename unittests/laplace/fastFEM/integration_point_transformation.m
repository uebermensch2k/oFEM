function point=integration_point_transformation(point)
%Affine transformation from the reference tetrahedra given by vertices
%[-1 -1 -1], [1 -1 -1], [-1 1 -1], [-1 -1 1]
% to the reference tetrahedra given by vertices (book of Solin on hp elem.)
%[0 0 0], [1 0 0], [0 1 0], [0 0 1]
switch size(point,2)
    case 3,
        point=([1 -1 -1; -1 1 -1; -1 -1 1]+1)\(point'+1);
    case 2,
        point=([1 -1; -1 1]+1)\(point'+1);
    otherwise, error('Only tringle transformation in 2D and tetrahedral transformation in 3D are implement.');
end




