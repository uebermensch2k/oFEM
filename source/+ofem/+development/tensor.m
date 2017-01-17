

classdef tensor
%     properties(Access=protected)
    properties(Access=public)
        A;
        N;
    end

    methods
        %%
        function obj=tensor(A)
        %TENSOR constructs tensor from a matrix
        %
            obj.A=A;
            obj.N=size(A,3);
        end

        %%
        function D=size(A,dim)
            D=size(A.A);
            if ~isempty(dim)
                D=D(dim);
            end
        end

        %%
        function A=ctranspose(A)
            A.A=conj(permute(A.A,[2,1,3]));
        end

        %%
        function A=transpose(A)
            A.A=permute(A.A,[2,1,3]);
        end

        %%
        function C=plus(A,B)
            C=ofem.tensor(A.A+B.A);
        end

        %%
        function C=minus(A,B)
            C=ofem.tensor(A.A-B.A);
        end

        %%
        function C=uplus(A)
            C=ofem.tensor(+A.A);
        end

        %%
        function C=uminus(A)
            C=ofem.tensor(-A.A);
        end

        %%
        function C=times(A,B)
            C=ofem.tensor(A.A.*B.A);
        end

        %%
        function C=mtimes(A,B)
            if size(A,3)~=size(B,3)
                error('MATLAB:ofem:tensor:mtimes:InvalidArgument',...
                      'Tensors must have the same size');
            end

            if isa(B,'ofem.tensor')
                NA=size(A.A,1);
                NB=size(B.A,2);

                B=B';

                C=ofem.tensor(reshape(dot(repmat (A.A,NB,1,1),...
                                          repelem(B.A,NA,1,1),...
                                          2),...
                                      NA,NB,[])...
                         );
            elseif isnumeric(B)
                C=A*ofem.tensor(repmat(B,1,1,size(A.A,3)));
            else
                error('MATLAB:ofem:tensor:mtimes:InvalidArgument',...
                      'Second argument must either be a tensor of a numeric array');
            end
        end

        %%
        function C=rdivide(A,B)
            C=ofem.tensor(A.A./B.A);
        end

        %%
        function C=ldivide(A,B)
            C=ofem.tensor(A.A.\B.A);
        end

        %%
        function C=power(A,B)
            C=ofem.tensor(A.A.^B.A);
        end

        %%
        function C=lt(A,B)
            C=ofem.tensor(A.A<B.A);
        end

        %%
        function C=gt(A,B)
            C=ofem.tensor(A.A>B.A);
        end

        %%
        function C=le(A,B)
            C=ofem.tensor(A.A<=B.A);
        end

        %%
        function C=ge(A,B)
            C=ofem.tensor(A.A>=B.A);
        end

        %%
        function C=ne(A,B)
            C=ofem.tensor(A.A~=B.A);
        end

        %%
        function C=eq(A,B)
            C=ofem.tensor(A.A==B.A);
        end

        %%
        function C=and(A,B)
            C=ofem.tensor(A.A&B.A);
        end

        %%
        function C=or(A,B)
            C=ofem.tensor(A.A|B.A);
        end

        %%
        function C=not(A)
            C=ofem.tensor(~A.A);
        end

        %%
        function disp(A)
            display(A.A);
        end

        %%
        function C=horzcat(A,B)
            C=ofem.tensor([A.A, B.A]);
        end

        %%
        function C=vertcat(A,B)
            C=ofem.tensor([A.A; B.A]);
        end

        %%
        function A=subsasgn(A,s,B)
            A.A=subsasgn(A.A,s,B);
        end

        %%
        function C=subsref(A,s)
            C=ofem.tensor(subsref(A.A,s));
        end

    end
end