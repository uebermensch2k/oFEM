

classdef elastic < handle
    properties(Access=protected)
        mesh;
        felem;
        lambda;
        mu;
    end

    methods(Access=protected,Static)
        %%
        function M=mass(pipj,detD,el,co)
        %MASS returns the mass matrix.
        %
        % M=MASS(pipj,detD,el,co) returns the mass matrix for the local set
        % of elements specified through el. pipj and detD are assumed to be
        % computed beforehand.
        %
        % see also OFEM.MESH.JACOBIANDATA, OFEM.FINITEELEMENT.PHIIPHIJ
        %
            Ns = size(pipj,1);
            Ne = size(el,1);

            Nc = size(co,1);
            Nd = size(co,2);
            
            
            % pipj for every coordiante dimension
            pipj=repelem(pipj,Nd,1);
            
            %checkerboard structure of the local mass matrix
            P=zeros(Nd*Ns, Nd*Ns);
            for i=1:Nd
                P(i:Nd:Nd*Ns,i:Nd:Nd*Ns) = pipj(i:Nd:end,:);
            end
            %the correct determinants meet the correct values
            M = detD*P(:)';
            %and now checkerboard it again
            M = M';
            
            I = repmat(1:Ns,Nd*Nd*Ns,1    ); I=I(:);
            J = repmat(1:Ns,Nd      ,Nd*Ns); J=J(:);
 
            I = Nd*el(:,I)-repmat(kron((Nd-1):-1:0,ones(1,Nd*Ns)),Ne,Ns); I=I';
            J = Nd*el(:,J)-repmat(kron(ones(1,Nd*Ns),(Nd-1):-1:0),Ne,Ns); J=J';

            M=sparse(I(:),J(:),M(:),Nd*Nc,Nd*Nc);
        end
        
        %%
        function S=stiffness(lambda,mu,DinvT,detD,dphi,w,el,co)
        %STIFFNESS returns the stiffness matrix.
        %
        % S=STIFFNESS(DinvT,detD,dphi,w,el,co) returns the stiffness matrix
        % for the local set of elements specified through el. The
        % computation is carried out in terms of a Gaussian quadrature, the
        % gradients dphi of the shape functions are expected to be computed
        % at the quadrature points with its third dimension associated to
        % the distinct quadrature points and W the weights of the
        % quadrature.
        %
        % see also OFEM.MESH.JACOBIANDATA, OFEM.FINITEELEMENT.DPHI,
        % OFEM.FINITEELEMENT.QUADDATA
        %

            Ns = size(dphi,1);
            Nd = size(dphi,2);
            Nq = size(dphi,3);
            Ne = size(el,1);
            Nc = size(co,1);

            I = repmat(1:Nd*Ns,Nd*Ns,1); I=I(:);
            J = repmat(1:Nd*Ns,1,Nd*Ns); J=J(:);

            % lambda = (nu*E)/((1+nu)*(1-2*nu)); %for plane strain
            % lambda = (nu*E)/((1+nu)^2);        %for plane stress
            % mu     = E/(2*(1+nu))

            % Approximated values for stainless steel
            % E = 200000000000 [pa] in SI units!
            % nu = 0.3

%             mu     = 76.9231;
%             lambda = 115.3846; 
%             mu=1; lambda=1;
%             lambda = 1.1538e+11; mu = 7.6923e+10;

            if ~isscalar(lambda) || ~isscalar(mu)
                
            end

            s1 = sqrt(  lambda+2*mu);
            s2 = sqrt(  lambda+  mu);
            s3 = sqrt(3*lambda+2*mu);

            sigma = [s1, lambda/s1         , lambda/s1              ; ...
                     0   , 2*sqrt(mu)*s2/s1, lambda*sqrt(mu)/(s1*s2); ...
                     0   , 0               , sqrt(mu)*s3/s2         ];

            idx   = [1,2; ...
                     1,3; ...
                     2,3];

            S=zeros(Nd*Ns*Nd*Ns,1,Ne);
            R=zeros(Nd*Ns,Ns,Ne);
            for k=1:Nq
                % DinvT is Nd*NexNd matrix
                % dphi is NsxNdxNl matrix
                % dphi(:,:,k)*DinvT' is NsxNd*Ne
                globdphi=reshape(dphi(:,:,k)*DinvT',Ns,Nd,Ne);

                % R is Nd*NsxNsxNe
                for i=1:Nd
                    for j=i:Nd
                        R(j:Nd:end,i,:) = globdphi(:,j,:)*sigma(i,j);
                    end
                end
                for i=1:nchoosek(Nd,2)
                    R(idx(i,1):Nd:end,i+Nd,:) = globdphi(:,idx(i,2),:)*sqrt(mu);
                    R(idx(i,2):Nd:end,i+Nd,:) = globdphi(:,idx(i,1),:)*sqrt(mu);
                end

                % globdphi is now NsxNdxNe matrix => extend along first
                % dimension and perform a dot product along the second
                S=S+w(k)*dot(R(I,:,:),R(J,:,:),2);
            end

            detD = repmat(detD',Nd*Ns*Nd*Ns,1);
            S    = squeeze(S).*detD;

            I = repmat(1:Ns,Nd*Nd*Ns,1    ); I=I(:);
            J = repmat(1:Ns,Nd      ,Nd*Ns); J=J(:);

            I = Nd*el(:,I)-repmat(kron((Nd-1):-1:0,ones(1,Nd*Ns)),Ne,Ns); I=I';
            J = Nd*el(:,J)-repmat(kron(ones(1,Nd*Ns),(Nd-1):-1:0),Ne,Ns); J=J';

            S=sparse(I(:),J(:),S(:),Nd*Nc,Nd*Nc);
        end

        %%
        function b=load(detD,phi,w,l,f,el,co)
        %LOAD returns the load vector.
        %
        % b=LOAD(detD,phi,w,l,f,el,co) computes the load vector in terms of
        % a Gaussian quadrature. w and l are the weights and quadrature
        % points of the rule, respectively. phi are the shape functions
        % evaluated at the quadrature points and f a functions handle
        % returning the value of the right hand side at arbitrary points.
        %
        % see also OFEM.MESH.JACOBIANDATA, OFEM.FINITEELEMENT.PHI,
        % OFEM.FINITEELEMENT.QUADDATA
        %
            Ne       = size(el,1);
            Nd       = size(co,2);
            Nn       = size(el,2);
            Nq       = size(w, 1);
            Nl       = size(l,2);
            Ns       = size(phi,1);
            F        = zeros(Nd*Ns*Ne,1);

            elco = reshape(permute(reshape(co(el(:,1:Nl),:),Ne,Nl,Nd),[1,3,2]),Nd*Ne,Nl);

            for p=1:Nq
                P    = reshape(elco*l(p,:)', Ne, Nd);
                fP   = reshape(f(P),[],1);
                F    = F+w(p)*kron(phi(:,p),fP);
            end
            F    = F.*repmat(detD,Ns*Nd,1);
            I    = repelem(1:Nn,1,Nd); I=I(:);
            I    = Nd*el(:,I) - repmat((Nd-1):-1:0,Ne,Ns); %I=I';
            b    = accumarray(I(:),F);
        end

        function b=neumann(meas,phi,w,l,g,faces,normals,co)
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
            Nl     = size(l  ,2);
            Nd     = size(co ,2);
            Nq     = size(w  ,1);
            Ns     = size(phi,1);
            Nf     = size(faces,1);
            Nn     = size(faces,2);

            F      = zeros(Nd*Ns*Nf,1);

            % faceco*l gives global quadrature points => evaluate there
            faceco  = reshape(permute(reshape(co(faces(:,1:Nl),:),Nf,Nl,Nd),[1,3,2]),Nd*Nf,Nl);
            normals = repelem(normals,Nd,1);
%             normals

            for q=1:Nq
                X    = reshape(faceco*l(q,:)',Nf,Nd);
                F    = F+w(q)*kron(phi(:,q),reshape(g(X,normals),[],1));
            end

            F    = F.*repmat(meas,Ns*Nd,1);
            I    = repelem(1:Nn,1,Nd); I=I(:);
            I    = Nd*faces(:,I) - repmat((Nd-1):-1:0,Nf,Ns); %I=I';
            b    = sparse(I(:),1,F,Nd*size(co,1),1);
%             b    = accumarray(faces(:),F);
        end
    end

    methods
        %%
        function obj=elastic(mesh,felem)
%             obj@ofem.laplace(mesh,felem);
            obj.mesh  = mesh;
            obj.felem = felem;
        end

        function setmaterial(obj,lambda,mu)
            obj.lambda=lambda;
            obj.mu=mu;
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
        % Note that the DOFs are dim x Nco. The solution is therefore
        % orderer as [u1T1;u2T1;u3T1;u1T2;u2T2;u3T2;...]. 
        %

            if nargin==1
                opt.S=1;
                opt.M=1;
                opt.b=0;
                opt.N=0;
            else
                opt=varargin{1};
                if ~isstruct(opt)
                    error('opt is expected to be a structure.');
                else

                if ~isfield(opt,'M'); opt.M=0; end
                if ~isfield(opt,'S'); opt.S=0; end
                if ~isfield(opt,'b'); opt.b=0; end
                if ~isfield(opt,'N'); opt.N=0; end
                end
            end

%             if opt.M==1
%                 warning('No mass matrix computation supported so far!');
%                 opt.M=0; % no stiffness matrix so far!
%             end

            if opt.b==1
                if ~isfield(opt,'load')
                    error('Load vector requested, but no function handle given!');
                elseif ~isa(opt.load,'function_handle')
                    error('load is expected to be a function handle!');
                end
            end
            
            if opt.N==1
                if ~isfield(opt,'neumann')
                    error('Neumann boundary requested, but no function handle given!');
                elseif ~isa(opt.neumann,'function_handle')
                    error('neumann is expected to be a function handle!');
                end

                if ~isfield(opt,'neumannidx')
                    error('You need to specifiy the Neumann boundary index!');
                end
            end

            tic;
            [DinvT,detD] = obj.mesh.jacobiandata;
            Np           = size(obj.mesh.parts,2);

            %% mass
            if opt.M==1
                M    = cell(Np,1);
                pipj = obj.felem.phiiphij(obj.mesh.dim);
            end

            if opt.S==1 || opt.b==1
               [w,l] = obj.felem.quaddata(obj.mesh.dim);
            end

            %% stiffness
            if opt.S==1
                S    = cell(Np,1);
                dphi = obj.felem.dphi(l);
            end
            
            b = sparse(obj.mesh.dim*size(obj.mesh.co,1),1);

            %% load vector
            if opt.b==1
                phi  = obj.felem.phi(l);
                b    = obj.load(detD,phi,w,l,opt.load,obj.mesh.el,obj.mesh.co);
            end
            
            %% Neumann boundary
            if opt.N==1
                [wN,lN] = obj.felem.quaddata(obj.mesh.dim-1);
                phiN   = obj.felem.phi(lN);
                [meas,faces,normals] = obj.mesh.neumann(opt.neumannidx);
                meas    = vertcat(meas {:});
                faces   = vertcat(faces{:});
                normals = vertcat(normals{:});
                b = b + obj.neumann(meas,phiN,wN,lN,opt.neumann,faces,normals,obj.mesh.co);
            end
            
            %Stiffnes and mass matrizes over the element parts
            if opt.S==1 || opt.M==1
                for i =1:Np
                    pIdx     = obj.mesh.parts{3,i};
                    elemsLoc = obj.mesh.el(pIdx,:);
                    detDLoc  = detD(pIdx);
                    
                    if opt.S==1
                        Didx     = false(1,length(pIdx)); Didx(pIdx)=true;
                        Didx     = repmat(Didx,obj.mesh.dim,1);
                        DinvTLoc = DinvT(Didx,:);
                        S{i}     = obj.stiffness(obj.lambda(i), obj.mu(i),DinvTLoc,detDLoc,dphi,w,elemsLoc,obj.mesh.co);
                    end
                    if opt.M==1
                        M{i}     = obj.mass(pipj,detDLoc,elemsLoc,obj.mesh.co);
                    end
                end
            end
            %% auxiliary variable
            info.time2assemble = toc;

            %% computed structure
            if opt.S==1; asm.S=S; end
            if opt.M==1; asm.M=M; end
            if opt.b==1 || opt.N==1; asm.b=b; end

            aux.detD=detD;
        end
        
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
                        n=(d+3)*(d+4)/2; % degree of point considering for least-squares
                        
                        % statistics and machine learning toolbox
                        co_idx = knnsearch(obj.mesh.co,obj.mesh.co,'K',n)';
                        
                        ui = u          (co_idx(:)  );
                        xi = obj.mesh.co(co_idx(:),1);
                        yi = obj.mesh.co(co_idx(:),2);
                        
                        clear co_idx;
                        
                        Ndof = size(obj.mesh.co,1);
                        
                        A = [ones(Ndof*n,1), xi, yi, xi.*yi, xi.^2, yi.^2];
                        
                        ui_rec = cellfun(@(A,b) A\b      , ...
                            mat2cell(A ,n*ones(1,Ndof),m), ...
                            mat2cell(ui,n*ones(1,Ndof),1), ...
                            'UniformOutput',false);
                        ui_rec = reshape(cell2mat(ui_rec),m,[])';
                        
                        clear A;
                        
                        x = obj.mesh.co(:,1);
                        y = obj.mesh.co(:,2);
                        
                        grad = [ ui_rec(:,2)+ui_rec(:,4).*y+2*ui_rec(:,5).*x, ...
                            ui_rec(:,3)+ui_rec(:,4).*x+2*ui_rec(:,6).*y ];
                        
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
                        
                        ui_rec = cellfun(@(A,b) A\b      , ...
                            mat2cell(A ,n*ones(1,Ndof),m), ...
                            mat2cell(ui,n*ones(1,Ndof),1), ...
                            'UniformOutput',false);
                        ui_rec = reshape(cell2mat(ui_rec),m,[])';
                        
                        clear A;
                        
                        x = obj.mesh.co(:,1);
                        y = obj.mesh.co(:,2);
                        z = obj.mesh.co(:,3);
                        
                        grad = [ ui_rec(:,2)+ui_rec(:,5).*y+ui_rec(:,6).*z+2*ui_rec(:,8).*x, ...
                            ui_rec(:,3)+ui_rec(:,5).*x+ui_rec(:,7).*z+2*ui_rec(:,9).*y,...
                            ui_rec(:,4)+ui_rec(:,6).*x+ui_rec(:,7).*y+2*ui_rec(:,10).*z];
                        
                        
                        % warning('3D gradient recovery not implemented yet!');
                end
            otherwise
                error('MATLAB:ofem:laplace:gradu:NotSupported',...
                          'Reconstruction of P2 gradient not implemented, yet!');
        end
        end
        
        %% Deformation Gradient, Strains and Stresses
        function F=defGrad(obj,u)
            % DEFGRAD computes the deformation gradient tensor at the DOFs
            %
            % The tensor is computed by calculating the coordiante
            % derivatives of the deformation and then adding the unit
            % tensor. All is performed in Voigt notation: F = I + DU/DX
            %
            
            F = zeros(obj.mesh.dim*size(obj.mesh.co,1),obj.mesh.dim);
            for i=1:obj.mesh.dim
                F(i:obj.mesh.dim:end,:) = obj.gradu(u(:,i));
            end
            for i=1:obj.mesh.dim
                F(i:obj.mesh.dim:end,i) = F(i:obj.mesh.dim:end,i)+1;
            end
        end
        
        function [E,e,s] = StrainStress(obj,u)
            
            co     = obj.mesh.co;
            dim    = obj.mesh.dim;
            F      = obj.defGrad(u);
            Np     = size(obj.mesh.parts,2);
            s      = cell(Np,1);
            e      = cell(Np,1);
            
            %% Green's Stain Tensor computation ( E =1/2*(F^T F - I))
            K= permute(reshape(F,dim,size(co,1),dim),[1 3 2]);
            I=repmat(1:dim,dim,1);
            J=I';
            EE=reshape(dot(K(I,:,:), K(J,:,:),2),dim,dim,[]);
            %this is the strain now, but in the wrong format so we have to change this
            E=0.5*(reshape(permute(EE,[1,3,2]),[],dim)-repmat(eye(dim,dim),...
                size(co,1),1));
            K= permute(reshape(E,dim,size(co,1),dim),[1 3 2]);
            

            E = reshape(K,dim*dim, size(u,1))';
            switch dim
                case 3
                     %isotropic stress tensor: only the idependen components need to be
                     %computed, this happens by using only the column components
                     %[1 5 9 6 3 2] (i.e. e_{11} e_{22} e_{33} e_{23} e_{13} e_{12}) of every row
                    for i =1:Np
                        pIdx     = obj.mesh.parts{3,i};
                        elemsLoc = obj.mesh.el(pIdx,:);
                        nIdx     = unique(elemsLoc);
                        uPart    = u(nIdx, :);
                        e{i} = repmat([ 1, 1, 1, 2, 2, 2],size(uPart,1), 1).*...
                            E(nIdx, [1,5,9,6,3,2]);
                        C = [2*obj.mu(i) + obj.lambda(i), obj.lambda(i), obj.lambda(i), 0,  0,  0;
                            obj.lambda(i),2*obj.mu(i)+obj.lambda(i), obj.lambda(i), 0,  0,  0;
                            obj.lambda(i),obj.lambda(i),2*obj.mu(i)+obj.lambda(i), 0,  0,  0;
                            0,               0,           0,         obj.mu(i),  0,  0;
                            0,               0,           0,          0, obj.mu(i),  0;
                            0,               0,           0,          0,  0, obj.mu(i)];
                        
                        s{i}=C*e{i}';
                    end
                case 2
                    warning('MATLAB:ofem:elastic:StrainStress:WordsOfWisdom',...
                            'You are calculating stresses in 2D, only plane stress is implemented so far, I will compute that')
                    %isotropic stress tensor: only the idependen components need to be
                    %computed, this happens by using only the column components
                    %[1 4 2] (i.e. e_{11} e_{22} e_{12}) of every row
                    for i =1:Np
                        pIdx     = obj.mesh.parts{3,i};
                        elemsLoc = obj.mesh.el(pIdx,:);
                        nIdx     = unique(elemsLoc);
                        uPart    = u(nIdx, :);
                        e{i} = repmat([ 1, 1,2],size(uPart,1), 1).*...
                            E(nIdx, [1,4,2]);
                        C = [2*obj.mu(i) + obj.lambda(i),  obj.lambda(i), 0;
                            obj.lambda(i),2*obj.mu(i)+obj.lambda(i),  0;
                             0,            0,                 obj.mu(i)];
                        
                        s{i}=C*e{i}';
                    end
                case 1
            end
            
        end
    end
end