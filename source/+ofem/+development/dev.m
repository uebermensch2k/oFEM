function d = dev(A, dim)
% DEV returns the deviator of a tensor
trA = ofem.tr(A,dim);
d=A-repmat(trA,1,dim*dim,1).*repmat(reshape(eye(dim,dim),1,[]),size(trA,1),1,size(A,3))/dim;
end
