

classdef matrixarray < double
%     properties(Access=protected)
%     properties(Access=public)
%     end

    methods
        %%
        function obj=matrixarray(data)
        %MATRIXARRAY constructs matrixarray from a matrix
        %
            if nargin==0
                data=0;
            end
            validateattributes(data,{'numeric','logical'},{'nonempty'});

            obj=obj@double(data);
        end

        %%
        function D=rot(A)
            narginchk(1,2);
%             validateattributes(A,{'ofem.matrixarray'},{});
            validateattributes(A(:,:,1),{'numeric'},{'vector','numel',2});

            dim=find(size(A(:,:,1))==2);

            if dim==1
                D=vertcat(-subsref(A,substruct('()',{2,':',':'})),...
                           subsref(A,substruct('()',{1,':',':'})));
            elseif dim==2
                D=horzcat(-subsref(A,substruct('()',{':',2,':'})),...
                           subsref(A,substruct('()',{':',1,':'})));
            else
                error('ofem:matrixarray:InvalidArgument',...
                      'dim can either be one or two');
            end
        end

        %%
        function C=dot(A,B,varargin)
            C=ofem.matrixarray(dot(double(A),double(B),varargin{:}));
        end

        %%
        function C=reshape(A,varargin)
            C=ofem.matrixarray(reshape(double(A),varargin{:}));
        end

        %%
        function C=cross(A,B,varargin)
            C=ofem.matrixarray(cross(double(A),double(B),varargin{:}));
        end

        %%
        function C=conj(obj)
            C=ofem.matrixarray(conj(double(obj)));
        end

        %%
        function C=ctranspose(obj)
        %CTRANSPOSE Point-wise complex transposition
            C=conj(permute(obj,[2,1,3]));
        end

        %%
        function C=transpose(obj)
        %CTRANSPOSE Point-wise transposition
            C=permute(obj,[2,1,3]);
        end

        %%
        function C=plus(A,B)
            %% only numeric or matrixarray supported
            narginchk(2,2);
            validateattributes(A,{'numeric'},{'nonempty'});
            validateattributes(B,{'numeric'},{'nonempty'});

            if isscalar(A) || isscalar(B) || (numel(size(A))==numel(size(B)) && all(size(A)==size(B)))
                C=ofem.matrixarray(double(A)+double(B));
            else
                if size(A,3)==size(B,3)
                    if isscalar(A(:,:,1))
                        C=ofem.matrixarray(repmat(double(A),size(B(:,:,1)))+double(B));
                    elseif isscalar(B(:,:,1))
                        C=ofem.matrixarray(double(A)+repmat(double(B),size(A(:,:,1))));
%                     % already asked for above
%                     elseif size(A(:,:,1))==size(B(:,:,1))
                    else
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                else % different lengths
                    if ~all(size(A(:,:,1))==size(B(:,:,1)))
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                    % then one must have length one
                    if size(A,3)==1
                        C=ofem.matrixarray(repmat(double(A),1,1,size(B,3))+double(B));
                    elseif size(B,3)==1
                        C=ofem.matrixarray(double(A)+repmat(double(B),1,1,size(A,3)));
                    else
                        error('ofem:matrixarray:LengthMismatch', ...
                              'Array must have same length or one must have length one');
                    end
                end
            end
        end

        %%
        function C=minus(A,B)
            C=plus(A,-B);
        end

        %%
        function C=uplus(A)
            C=ofem.matrixarray(+double(A));
        end

        %%
        function C=uminus(A)
            C=ofem.matrixarray(-double(A));
        end

        %%
        function C=times(A,B)
            %% only numeric or matrixarray supported
            narginchk(2,2);
            validateattributes(A,{'numeric'},{'nonempty'});
            validateattributes(B,{'numeric'},{'nonempty'});

            if isscalar(A) || isscalar(B) || (numel(size(A))==numel(size(B)) && all(size(A)==size(B)))
                C=ofem.matrixarray(double(A).*double(B));
            else
                if size(A,3)==size(B,3)
                    if isscalar(A(:,:,1))
                        C=ofem.matrixarray(repmat(double(A),size(B(:,:,1))).*double(B));
                    elseif isscalar(B(:,:,1))
                        C=ofem.matrixarray(double(A).*repmat(double(B),size(A(:,:,1))));
%                     % already asked for above
%                     elseif size(A(:,:,1))==size(B(:,:,1))
                    else
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                else % different lengths
                    if ~all(size(A(:,:,1))==size(B(:,:,1)))
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                    % then one must have length one
                    if size(A,3)==1
                        C=ofem.matrixarray(repmat(double(A),1,1,size(B,3)).*double(B));
                    elseif size(B,3)==1
                        C=ofem.matrixarray(double(A).*repmat(double(B),1,1,size(A,3)));
                    else
                        error('ofem:matrixarray:LengthMismatch', ...
                              'Array must have same length or one must have length one');
                    end
                end
            end
        end

        %%
        function C=mtimes(A,B)
            %% only numeric or matrixarray supported
            narginchk(2,2);
            validateattributes(A,{'numeric'},{'nonempty'});
            validateattributes(B,{'numeric'},{'nonempty'});

            %% one of the factors is scalar
            if isscalar(A)
                C=ofem.matrixarray(double(A)*double(B));
                return;
            elseif isscalar(B)
                C=ofem.matrixarray(double(A)*double(B));
                return;
            end

            %% one of the factors is scalar per entry
            if isscalar(A(:,:,1))
                NrB=size(B,1);
                NcB=size(B,2);
                C=repelem(A,NrB,NcB,1).*B;
                return;
            elseif isscalar(B(:,:,1))
                NrA=size(A,1);
                NcA=size(A,2);
                C=A.*repelem(B,NrA,NcA,1);
                return;
            end

            %% check for inner matrix dimension
            if size(A,2)~=size(B,1)
                error('ofem:matrixarray:InnerDimensionMismatch',...
                      'Inner matrix dimension must agree');
            end

            %% one of the factors is matrix
            if ismatrix(A)
                NrA=size(A,1);
                NrB=size(B,1);
                NcB=size(B,2);

                tmp=double(A)*reshape(double(B),NrB,[]);
                C=ofem.matrixarray(reshape(tmp,NrA,NcB,[]));
                return;
            elseif ismatrix(B)
                NrA=size(A,1);
                NcA=size(A,2);
                NcB=size(B,2);

                tmp=reshape(permute(double(A),[1,3,2]),[],NcA)*double(B);
                C=ofem.matrixarray(permute(reshape(tmp,NrA,[],NcB),[1,3,2]));
                return;
            end

            %% both factors are matrixarrays
            if size(A,3)~=size(B,3)
                error('ofem:matrixarray:InvalidArgument',...
                      'Arrays must have the same length');
            end

            NA=size(A,1);
            NB=size(B,2);

            B=B';

            tmp = dot(repmat(double(A),NB,1,1),repelem(double(B),NA,1,1),2);
            C=ofem.matrixarray(reshape(tmp,NA,NB,[]));
        end

        %%
        function C=rdivide(A,B)
            %% only numeric or matrixarray supported
            narginchk(2,2);
            validateattributes(A,{'numeric'},{'nonempty'});
            validateattributes(B,{'numeric'},{'nonempty'});

            if isscalar(A) || isscalar(B) || (numel(size(A))==numel(size(B)) && all(size(A)==size(B)))
                C=ofem.matrixarray(double(A)./double(B));
            else
                if size(A,3)==size(B,3)
                    if isscalar(A(:,:,1))
                        C=ofem.matrixarray(repmat(double(A),size(B(:,:,1)))./double(B));
                    elseif isscalar(B(:,:,1))
                        C=ofem.matrixarray(double(A)./repmat(double(B),size(A(:,:,1))));
%                     % already asked for above
%                     elseif size(A(:,:,1))==size(B(:,:,1))
                    else
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                else % different lengths
                    if ~all(size(A(:,:,1))==size(B(:,:,1)))
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                    % then one must have length one
                    if size(A,3)==1
                        C=ofem.matrixarray(repmat(double(A),1,1,size(B,3))./double(B));
                    elseif size(B,3)==1
                        C=ofem.matrixarray(double(A)./repmat(double(B),1,1,size(A,3)));
                    else
                        error('ofem:matrixarray:LengthMismatch', ...
                              'Array must have same length or one must have length one');
                    end
                end
            end
        end

        %%
        function C=ldivide(A,B)
            %% only numeric or matrixarray supported
            narginchk(2,2);
            validateattributes(A,{'numeric'},{'nonempty'});
            validateattributes(B,{'numeric'},{'nonempty'});

            if isscalar(A) || isscalar(B) || (numel(size(A))==numel(size(B)) && all(size(A)==size(B)))
                C=ofem.matrixarray(double(A).\double(B));
            else
                if size(A,3)==size(B,3)
                    if isscalar(A(:,:,1))
                        C=ofem.matrixarray(repmat(double(A),size(B(:,:,1))).\double(B));
                    elseif isscalar(B(:,:,1))
                        C=ofem.matrixarray(double(A).\repmat(double(B),size(A(:,:,1))));
%                     % already asked for above
%                     elseif size(A(:,:,1))==size(B(:,:,1))
                    else
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                else % different lengths
                    if ~all(size(A(:,:,1))==size(B(:,:,1)))
                        error('ofem:matrixarray:DimensionMismatch', ...
                              'Each array entry must either be a scalar or dimensions must match');
                    end
                    % then one must have length one
                    if size(A,3)==1
                        C=ofem.matrixarray(repmat(double(A),1,1,size(B,3)).\double(B));
                    elseif size(B,3)==1
                        C=ofem.matrixarray(double(A).\repmat(double(B),1,1,size(A,3)));
                    else
                        error('ofem:matrixarray:LengthMismatch', ...
                              'Array must have same length or one must have length one');
                    end
                end
            end
        end


%         %%
%         function C=power(A,B)
%             %% only numeric or matrixarray supported
%             narginchk(2,2);
%             validateattributes(A,{'numeric','ofem.matrixarray'},{'3d','nonempty'});
%             validateattributes(B,{'numeric','ofem.matrixarray'},{'3d','nonempty'});
% 
%             %% one of the addends is scalar
%             if isscalar(A)
%                 C=ofem.matrixarray(A.^B.A);
%                 return;
%             elseif isscalar(B)
%                 C=ofem.matrixarray(A.A.^B);
%                 return;
%             end
% 
%             %% check for matrix dimension
%             if size(A,1)~=size(B,1) || size(A,2)~=size(B,2)
%                 error('MATLAB:ofem:matrixarray:power:InvalidArgument',...
%                       'Matrix dimensions must agree');
%             end
% 
%             %% one of the addends is matrix
%             if ismatrix(A)
%                 for i=1:size(A,1)
%                     for j=1:size(A,2)
%                         B.A(i,j,:)=B.A(i,j,:).^A(i,j);
%                     end
%                 end
%                 C=B;
%                 return;
%             elseif ismatrix(B)
%                 for i=1:size(A,1)
%                     for j=1:size(A,2)
%                         A.A(i,j,:)=A.A(i,j,:).^B(i,j);
%                     end
%                 end
%                 C=A;
%                 return;
%             end
% 
%             %% both factors are matrixarrays
%             if size(A,3)~=size(B,3)
%                 error('MATLAB:ofem:matrixarray:power:InvalidArgument',...
%                       'Arrays must have the same length');
%             end
% 
%             C=ofem.matrixarray(A.A.^B.A);
%         end
%
    end
end