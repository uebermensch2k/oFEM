function rhs = myload(P)
% my load returns the load vecto [f_x, f_y, f_z]'  as long as 3 times the
% coordinates

%    rhs = zeros(3*size(P,1),1);
%    rhs(3:3:end) = -0.003;
   rhs = ofem.matrixarray(zeros(size(P)));
   rhs(3,:,:) = -0.003;

end