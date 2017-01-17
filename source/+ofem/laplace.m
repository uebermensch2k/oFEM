

classdef laplace < handle
    properties(Access=protected)
        mesh;
        felem;
    end

    methods(Access=protected,Static)

        %%
        function M=mass(pipj,detD,el,co)
        %MASS returns the mass matrix.
        %
        % M=MASS(pipj,detD,el,co) returns the mass matrix for the local set
        % of elements specified through el.co are the coordinates of the
        % mesh. detD is the per element determinant of the Jacobian of the
        % diffeomorphism taking the reference element to the global one.
        % E.g., detD(1,1,i) it the determinant for the i-th element. See
        % OFEM.MESH.JACOBIANDATA on how these quantities are organized.
        % pipj is the integral of the reference element over
        % \phi_i\cdot\phi_j, with \phi_k being the k-th basis function. See
        % OFEM.FINITEELEMENT.PHIIPHIJ on how the data is orginized.
        %
        % See also OFEM.MESH.JACOBIANDATA, OFEM.FINITEELEMENT.PHIIPHIJ
        %
            Ns = size(pipj,1);
            Nc = size(co  ,3);

            M = pipj*detD;

            J = repmat(1:Ns,Ns,1);
            I = el(:,J )';
            J = el(:,J')';

            M=sparse(I(:),J(:),M(:),Nc,Nc);
        end

        %%
        function S=stiffness(DinvT,detD,dphi,w,el,co)
        %STIFFNESS returns the stiffness matrix.
        %
        % S=STIFFNESS(DinvT,detD,dphi,w,el,co) returns the stiffness matrix
        % for the local set of elements specified through el. co are the
        % coordinates of the mesh. DinvT and detD are, respectively, the
        % transposed inverse of the Jacobian of Phi and the Jacobians
        % determinant per element, where Phi denotes the diffeomorphism
        % taking the reference element to the global one. E.g., for the
        % i-th element DinvT(:,:,i) is the i-th Jacobian's transposed
        % inverse and detD(1,1,i) the determinant. See
        % OFEM.MESH.JACOBIANDATA on how these quantities are organized.
        %
        % The computation is carried out in terms of a Gaussian quadrature.
        % The gradients dphi at quadrature points and weights w of the
        % quadrature rule are expected as returned by
        % OFEM.FINITEELEMENT.DPHI. The quadrature points can be queried
        % from OFEM.FINITEELEMENT.QUADDATA.
        %
        % See also OFEM.MESH.JACOBIANDATA, OFEM.FINITEELEMENT.DPHI,
        % OFEM.FINITEELEMENT.QUADDATA
        %

            Ns = size(dphi,1);
            Nq = size(dphi,3);
            Ne = size(el,1);
            Nc = size(co,3);

%             tic
            S=ofem.matrixarray(zeros(Ns,Ns,Ne));
            for q=1:Nq
                globdphi=DinvT*dphi(:,:,q)';
                S=S+(globdphi'*globdphi)*w(q);
            end

            S=S*detD;

            J = repmat(1:Ns,Ns,1);
            I = el(:,J )';
            J = el(:,J')';

            S = sparse(I(:),J(:),S(:),Nc,Nc);
%             toc
        end

%         %%
%         function N=neumann_bd(DinvT,normals,phi,dphi,detL,ss,el,co)
%         % NEUMANN_BD
%         %
%         % DinvT is Nd*Ne x Nd matrix
%         % normals is Ne x Nd matrix
%         % piDpj is Ns x Ns x Nd array
%         % detD is Ne vector
%         %
%         
%             Ns = size(dphi,1);
%             Nd = size(dphi,2);
%             Nq = size(dphi,3);
%             Ne = size(el,1);
%             Nc = size(co,1);
% 
%             % since S is symmetric we need only lower left triangular part
%             idx = tril(true(Ns));
%             I   = repmat(1:Ns,Ns,1); J = I'    ;
%             I   = I(idx)           ; J = J(idx);
% 
%             S=zeros(Ns*(Ns+1)/2,1,Ne);
%             for q=1:Nq
%                 % DinvT is Nd*NexNd matrix
%                 % dphi is NsxNdxNl matrix
%                 % dphi(:,:,k)*DinvT' is NsxNd*Ne
%                 globdphi=reshape(dphi(:,:,q)*DinvT',Ns,Nd,Ne);
% 
%                 % globdphi is now NsxNdxNe matrix => extend along first
%                 % dimension and perform a dot product along the second
%                 % NOTE: Since the stiffness matrix ist symmetric, only
%                 % lower triangular part needs to be computed.
%                 S=S+w(q)*dot(globdphi(I,:,:),globdphi(J,:,:),2);
%             end
% 
%             S = squeeze(S).*repmat(detD',Ns*(Ns+1)/2,1);
% 
% 
%             I = el(:,I)';
%             J = el(:,J)';
% 
%             S = sparse(I(:),J(:),S(:),Nc,Nc);
%             S = S+S'-diag(diag(S));
% 
%             Nc    = size(co,1);
%             Nd    = size(co,2);
%             Nss   = numel(ss);
% 
%             N = cell(Nss,1);
%             
%             for k=1:Nss
%                 eI = ss{k};
%                 DI = false(1,length(pIdx)); DI(pIdx)=true; DI = repmat(DI,Nd,1); DI = DI(:);
%                 Ns = size(piDpj{k},1);
%                 Ne = numel(eI);
% 
%                 % multiply by detD
%                 tmp=normals(eI,:).*repmat(detD(eI,:),Nd);
%                 % compute n^T*D^-T => N is Nd x Ne vector
%                 tmp=reshape(dot(repelem(tmp,Nd,1),DinvT(DI,:),2),Nd,[],1);
%                 % extend => N is Ns*Ne x Ns x Nd array
%                 tmp=repelem(repmute(tmp,[2,3,1]),Ns,Ns);
%                 % compute \phi_i n^TD^-T \nabla \phi_j => N is Ns*Ne x Ne
%                 N{k}=dot(tmp,repmat(piDpj{k},Ne,1),3);
%             end
% 
%             N     = horzcat(N{:});
%             eIall = vertcat(ss{:});
% 
%             J = repmat(1:Ns,Ns,1); I = J'  ;
%             I = I(:)             ; J = J(:); 
% 
%             I = el(eIall,I)';
%             J = el(eIall,J)';
% 
%             N = sparse(I(:),J(:),N(:),Nc,Nc);
%         end

        %%
        function b=load(pIdx,detD,phi,w,l,f,el,co)
        %LOAD returns the load vector.
        %
        % b=LOAD(pIdx,detD,phi,w,l,f,el,co) computes the load vector in
        % terms of a Gaussian quadrature. w and l are the weights and
        % quadrature points of the rule, respectively. phi are the shape
        % functions evaluated at the quadrature points and f a functions
        % handle returning the value of the right hand side at arbitrary
        % points. pIdx is the index of the part of the mesh.
        %
        % see also OFEM.MESH.JACOBIANDATA, OFEM.FINITEELEMENT.PHI,
        % OFEM.FINITEELEMENT.QUADDATA
        %
            Nl   = size(l  ,2); % number of barycentric coordinates
            Nc   = size(co ,3);
            Nq   = size(w  ,1);
            Ns   = size(phi,1);
            Ne   = size(el ,1);

            elco = reshape(co(:,:,el(:,1:Nl)'),[],Nl,Ne);
            F    = ofem.matrixarray(zeros(1,Ns,Ne));
            for q=1:Nq
                X = elco*(l(q,:)');
                F = F + f(X,pIdx)*(w(q)*phi(:,q)');
            end

            F = F*detD;
            b = sparse(el(:),1,F(:),Nc,1);
        end

        function b=neumann(meas,phi,w,l,g,faces,co)
        %NEUMANN returns the Neumann vector.
        %
        % b=NEUMANN(meas,phi,w,l,g,faces,co) computes the Neumann vector in
        % terms of a Gaussian quadrature. w and l are the weights and
        % quadrature points of the rule, respectively. phi are the shape
        % functions evaluated at the quadrature points and g a functions
        % handle returning the Neumann data at arbitrary points.
        %
        % see also OFEM.MESH.JACOBIANDATA, OFEM.MESH.NEUMANN,
        % OFEM.FINITEELEMENT.PHI, OFEM.FINITEELEMENT.QUADDATA
        %
            
            Nl     = size(l  ,2); % number of barycentric coordinates
            Nc     = size(co ,3);
            Nq     = size(w  ,1);
            Ns     = size(phi,1);
            Nf     = size(faces,1);

            faceco = reshape(co(:,:,faces(:,1:Nl)'),[],Nl,Nf);
            F      = ofem.matrixarray(zeros(1,Ns,Nf));
            for q=1:Nq
                X = faceco*(l(q,:)');
                F = F + g(X)*(w(q)*phi(:,q)');
            end

            F = F*meas;
            b = sparse(faces(:),1,F(:),Nc,1);
        end
    end

    methods
        %%
        function obj=laplace(mesh,felem)
        %LAPLACE construct the object from a OFEM.MESH and a
        %OFEM.FINITEELEMENT class.
        %
        % see also OFEM.MESH, OFEM.FINITEELEMENT
        %
            obj.mesh=mesh;
            obj.felem=felem;
        end

        %%
        function [asm,info,aux] = assemble(obj,varargin)
        %ASSEMBLE assembles the desired matrices and load vector.
        %
        % [asm,info,aux]=assemble computes the stiffness and the mass
        % matrices for each part specified in the OFEM.MESH class
        % associated with this class.
        %
        % [asm,info,aux]=assemble(opt), additionally, makes it possible to
        % specify which matrix and/or load vector to compute. E.g. if
        % opt.S==1, opt.M==0, opt.b==0, then only the stiffness matrix is
        % computed.
        %

            if nargin==1
                % Laplace equation with homogenous Neumann boundary
                opt.S=1;
                opt.M=0;
                opt.b=0;
                opt.N=0;
                opt.D=0;
            else
                if nargin>3
                    warning('ofem:laplace:TooManyArguments',...
                            'I''m expecting at most two arguments. I''ll skip the rest of them!');
                end
                opt=varargin{1};

                if ~isstruct(opt)
                    error('ofem:laplace:InvalidArgument',...
                          'opt is expected to be a structure.');
                end

                if ~isfield(opt,'M'); opt.M=0; end
                if ~isfield(opt,'S'); opt.S=0; end
                if ~isfield(opt,'b'); opt.b=0; end
                if ~isfield(opt,'N'); opt.N=0; end
                if ~isfield(opt,'D'); opt.D=0; end
            end

            tic;
            Np    = size(obj.mesh.parts,2);
            [w,l] = obj.felem.quaddata(obj.mesh.dim);

            %% handle load vector
            if opt.b==1
                if ~isfield(opt,'load')
                    error('ofem:laplace:InvalidArgument',...
                          'Load vector requested, but no function handle given!');
                elseif ~isa(opt.load,'function_handle')
                    error('ofem:laplace:InvalidArgument',...
                          'load is expected to be a function handle!');
                end

                asm.b_l = cell(Np,1);
                phi     = obj.felem.phi(l);
            end

            %% handle Neumann boundary
            if opt.N==1
                if ~isfield(opt,'neumann')
                    error('ofem:laplace:InvalidArgument',...
                          'Neumann boundary requested, but no function handle given!');
                elseif ~isa(opt.neumann,'function_handle')
                    error('ofem:laplace:InvalidArgument',...
                          'neumann is expected to be a function handle!');
                end

                if ~isfield(opt,'neumannidx')
                    error('ofem:laplace:InvalidArgument',...
                          'You need to specifiy the Neumann boundary index!');
                end

                [wN,lN] = obj.felem.quaddata(obj.mesh.dim-1);
                phiN    = obj.felem.phi(lN);
                [meas,faces,~] = obj.mesh.neumann(opt.neumannidx);
                meas  = vertcat(meas {:});
                faces = vertcat(faces{:});
                asm.b_N = obj.neumann(meas,phiN,wN,lN,opt.neumann,faces,obj.mesh.co);
            end

            %% handle Dirichlet boundary
            if opt.D==1
                if ~isfield(opt,'dirichlet')
                    error('ofem:laplace:InvalidArgument',...
                          'Dirichlet boundary requested, but no function handle given!');
                elseif ~isa(opt.dirichlet,'function_handle')
                    error('ofem:laplace:InvalidArgument',...
                          'dirichlet is expected to be a function handle!');
                end

                if ~isfield(opt,'dirichletidx')
                    error('ofem:laplace:InvalidArgument',...
                          'You need to specifiy the Dirichlet boundary index!');
                end

                asm.b_D=[];
            end

            [DinvT,detD] = obj.mesh.jacobiandata;

            %% Handle mass matrix
            if opt.M==1
                asm.M = cell(Np,1);
                pipj  = obj.felem.phiiphij(obj.mesh.dim);
            end

            %% Handle stiffness matrix
            if opt.S==1
                asm.S = cell(Np,1);
                dphi  = obj.felem.dphi(l);
            end

            %% compute Stiffnes and mass matrizes and load vector over the element parts
            if opt.S==1 || opt.M==1 || opt.b==1
                for i=1:Np
                    pIdx     = obj.mesh.parts{3,i};
                    elemsLoc = obj.mesh.el(pIdx,:);
                    detDLoc  = detD(:,:,pIdx);

                    if opt.b==1
                        asm.b_l{i} = obj.load(pIdx,detDLoc,phi,w,l,opt.load,elemsLoc,obj.mesh.co);
                    end
                    if opt.S==1
                        asm.S{i}   = obj.stiffness(DinvT(:,:,pIdx),detDLoc,dphi,w,elemsLoc,obj.mesh.co);
                    end
                    if opt.M==1
                        asm.M{i}   = obj.mass(pipj,detDLoc,elemsLoc,obj.mesh.co);
                    end
                end
            end

            %% auxiliary variable
            info.time2assemble = toc;

%             %% computed structure
%             if opt.S==1; asm.S=S; end
%             if opt.M==1; asm.M=M; end
%             if opt.b==1 || opt.N==1; asm.b=b; end
% 
%             aux.detD=detD;
            aux=[];
        end

        %%
        function grad=gradu(obj,u)
        %GRADU computes the gradient at DOFs.
        %
        % grad=gradu(u) computes the gradient grad of the FEM solution u at
        % DOFs. grad is a Ndofs by Nd matrix, where Ndofs are the number of
        % nodes and Nd the dimension of the spatial space.
        %
        
        switch obj.felem
            case ofem.finiteelement.P1
                
                switch obj.mesh.dim
                    case 2
                        % d denotes the dimension of polynomial space,
                        % i.e. for P1 elements it is 1
                        d=1;             % degree of finite element space
                        m=(d+2)*(d+3)/2; % polynomial degree of approximant

                        Ndof = size(obj.mesh.co,3);
                        
                        % statistics and machine learning toolbox
                        co = double(squeeze(obj.mesh.co))';
                        I  = obj.mesh.el(:,[1 1 2 2 3 3]);
                        J  = obj.mesh.el(:,[2 3 1 3 1 2]);
                        h  = co(I,:)-co(J,:); h=sqrt(dot(h,h,2));
                        h  = accumarray(I(:),h(:),[],@max,[],true);
                        ht = full(h);
                        D  = pdist2(co,co);
                        while 1
                            idx  = D<repmat(ht,1,Ndof);
                            idx2 = sum(idx,2)<m;
                            if ~any(idx2)
                                break;
                            end
                            ht(idx2)=ht(idx2)+h(idx2);
                        end
                        clear I J idx2 maxR maxRt D;

                        ui_r=zeros(6,1,Ndof);
                        for i=1:Ndof
                            idxl=find(idx(i,:));
                            ui = u (idxl  );
                            xi = co(idxl,1);
                            yi = co(idxl,2);

                            A  = [ones(numel(idxl),1), xi, yi, xi.*yi, xi.^2, yi.^2];

                            ui_r(:,:,i) = A\ui;
                        end

%                         n=(d+3)*(d+4)/2; % degree of point considering for least-squares
%                         co_idx = knnsearch(co,co,'K',n)';
% 
%                         ui = reshape(u (co_idx(:)  ),n,1,[]);
%                         xi = reshape(co(co_idx(:),1),n,1,[]);
%                         yi = reshape(co(co_idx(:),2),n,1,[]);
%                         
%                         clear co_idx;
% 
%                         A = [ones(n,1,Ndof), xi, yi, xi.*yi, xi.^2, yi.^2];
% 
%                         ui_r=zeros(6,1,Ndof);
%                         for i=1:Ndof
%                             ui_r(:,:,i) = A(:,:,i)\ui(:,:,i);
%                         end

                        x = obj.mesh.co(1,:,:);
                        y = obj.mesh.co(2,:,:);

                        grad = [ ui_r(2,:,:)+ui_r(4,:,:).*y+2*ui_r(5,:,:).*x, ...
                                 ui_r(3,:,:)+ui_r(4,:,:).*x+2*ui_r(6,:,:).*y ];

                        grad = reshape(grad,2,[])';

%                         ui_r = cellfun(@(A,b) A\b      , ...
%                             mat2cell(A ,n*ones(1,Ndof),m), ...
%                             mat2cell(ui,n*ones(1,Ndof),1), ...
%                             'UniformOutput',false);
%                         ui_r = reshape(cell2mat(ui_r),m,[])';
%                         
%                         clear A;
%                         
%                         x = co(:,1);
%                         y = co(:,2);
%                         
%                         grad = [ ui_r(:,2)+ui_r(:,4).*y+2*ui_r(:,5).*x, ...
%                                  ui_r(:,3)+ui_r(:,4).*x+2*ui_r(:,6).*y ];
                        
                    case 3
                        % d denotes the dimension of polynomial space,
                        % i.e. for P1 elements it is 1
                        %d=1;             % degree of finite element space
                        %m=(d+2)*(d+3)/2; % polynomial degree of approximant
                        m = 10;
                        % n=(d+3)*(d+4)/2; % degree of point considering for least-squares
                        n = 3*m;
                        
                        % statistics and machine learning toolbox
                        co_idx = knnsearch(obj.mesh.co,obj.mesh.co,'K',n)';
                        
                        %                 ui = u (co_idx(:)  );
                        %                 xi = co(co_idx(:),1);
                        %                 yi = co(co_idx(:),2);
                        %                 zi = co(co_idx(:),3);
                        
                        ui = u          (co_idx(:)  );
                        xi = obj.mesh.co(co_idx(:),1);
                        yi = obj.mesh.co(co_idx(:),2);
                        zi = obj.mesh.co(co_idx(:),3);
                        
                        
                        
                        clear co_idx;
                        
                        Ndof = size(obj.mesh.co,1);
                        
                        A = [ones(Ndof*n,1), xi, yi, zi, xi.*yi ,xi.*zi, yi.*zi, xi.^2, yi.^2, zi.^2];
                        
                        ui_r = cellfun(@(A,b) A\b      , ...
                            mat2cell(A ,n*ones(1,Ndof),m), ...
                            mat2cell(ui,n*ones(1,Ndof),1), ...
                            'UniformOutput',false);
                        ui_r = reshape(cell2mat(ui_r),m,[])';
                        
                        clear A;
                        
                        x = obj.mesh.co(:,1);
                        y = obj.mesh.co(:,2);
                        z = obj.mesh.co(:,3);
                        
                        grad = [ ui_r(:,2)+ui_r(:,5).*y+ui_r(:,6).*z+2*ui_r(:,8).*x, ...
                            ui_r(:,3)+ui_r(:,5).*x+ui_r(:,7).*z+2*ui_r(:,9).*y,...
                            ui_r(:,4)+ui_r(:,6).*x+ui_r(:,7).*y+2*ui_r(:,10).*z];
                        
                        
                        % warning('3D gradient recovery not implemented yet!');
                end
            otherwise
                error('ofem:laplace:NotSupported',...
                      'Reconstruction of P2 gradient not implemented, yet!');
        end
        end

        %%
        function plot(obj,u,varargin)
        %PLOT plots the solution.
        %
        co = double(squeeze(obj.mesh.co)');

        name = 'u';
        if nargin==3
            name=varargin{1};
        end

        switch obj.mesh.dim
            case 2
                trimesh(obj.mesh.el          , ...
                        co(:,1)     , ...
                        co(:,2)     , ...
                        u                    , ...
                        'FaceColor', 'interp', ...
                        'EdgeColor', 'none'  );
                xlabel('x-axis');
                ylabel('y-axis');
                zlabel(name);
                colorbar;

            case 3
                tetramesh(obj.mesh.el, ...
                          co         , ...
                          u          , ...
                          'FaceColor', 'interp', ...
                          'EdgeColor', 'none'  );
%                 trimesh(obj.mesh.el          , ...
%                         co(:,1)     , ...
%                         co(:,2)     , ...
%                         co(:,3)     , ...
%                         u                    , ...
%                         'FaceColor', 'interp', ...
%                         'EdgeColor', 'none'  );
                xlabel('x-axis');
                ylabel('y-axis');
                zlabel(name);
                colorbar;
        end
        end
    end
end