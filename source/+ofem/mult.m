function C=mult(A,B)
%MULT retruns the matrixproduct of two tensors written in
%           Voigt noation, i.e., if a 2nd order tensor is
%           identified with a matrix it returns the matrix
%           product. Tensormult operates on all element,
%           meaning it retruns the product over the 3rd
%           dimension
sA=size(A);
sB=size(B);

if sA(3)==1 && sB(3)==1 % just matrices
    C=A*B;
elseif isscalar(A) || isscalar(B)
    C=A*B;
% elseif sA(3)==1 || sB(3)==1
elseif numel(sA)<=2 || numel(sB)<=2
    if sA(2)~=sB(1)
        error('MATLAB:ofem:tensormult:InvalidArgument',...
            'Inner dimension must agree');
    end
    
%     if sA(3)==1
    if numel(sA)<=2
        C=ofem.mult(repmat(A,1,1,sB(3)),B);
    else
        C=ofem.mult(A,repmat(B,1,1,sA(3)));
    end
elseif sA(3)==sB(3)
    if prod(sA(1:2))==1
        C=repmat(A,sB(1),sB(2),1).*B;
    elseif prod(sB(1:2))==1
        C=A.*repmat(B,sA(1),sA(2),1);
    else
        NA=size(A,1);
        NB=size(B,2);
        C=reshape(dot(repmat(A,NB,1,1),repelem(permute(B,[2,1,3]),NA,1,1),2),NA,NB,[]);
    end
else
    error('MATLAB:ofem:tensormult:InvalidArgument',...
        'size(A,3)~=size(B,3)');
end
end
