

classdef elastic < handle
    properties(Access=protected)
        mesh;
        fe;
        qr;
        lambda;
        mu;
    end

    methods(Access=protected,Static)
        %%
        function M=mass(c, detD,pipj,el,co)
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

            Nc = size(co,3);
            Nd = size(co,1);
            
            
            % pipj for every coordiante dimension
            pipj=repelem(pipj,Nd,1);
            
            %checkerboard structure of the local mass matrix
            P=zeros(Nd*Ns, Nd*Ns);
            for i=1:Nd
                P(i:Nd:Nd*Ns,i:Nd:Nd*Ns) = pipj(i:Nd:end,:);
            end
            %the correct determinants meet the correct values
            M = detD*c*P(:)';
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
            Nc = size(co,3);

%             I = repmat(1:Nd*Ns,Nd*Ns,1); I=I(:);
%             J = repmat(1:Nd*Ns,1,Nd*Ns); J=J(:);

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
             
            S     = ofem.matrixarray(zeros(Ns*Nd, Ns*Nd, Ne));
            R     = ofem.matrixarray(zeros(Ns, Nd*Ns,Ne));
            for k=1:Nq
                % DinvT is Nd*NexNd matrix
                % dphi is NsxNdxNl matrix
                % dphi(:,:,k)*DinvT' is NsxNd*Ne
                globdphi = DinvT*dphi(:,:,k)';
                
                % R is Nd*NsxNsxNe
                for i=1:Nd
                    for j=i:Nd
                        R(i,j:Nd:end,:) = globdphi(j,:,:)*sigma(i,j);
                    end
                end
                for i=1:nchoosek(Nd,2)
                     R(i+Nd,idx(i,1):Nd:end,:) = globdphi(idx(i,2),:,:)*sqrt(mu);
                     R(i+Nd,idx(i,2):Nd:end,:) = globdphi(idx(i,1),:,:)*sqrt(mu);
                 end

                % globdphi is now NsxNdxNe matrix => extend along first
                % dimension and perform a dot product along the second
                S =S+R'*w(k)*R;
            end

            S=S*detD;

            I = repmat(1:Ns,Nd*Nd*Ns,1    ); I=I(:);
            J = repmat(1:Ns,Nd      ,Nd*Ns); J=J(:);

            I = Nd*el(:,I)-repmat(kron((Nd-1):-1:0,ones(1,Nd*Ns)),Ne,Ns); I=I';
            J = Nd*el(:,J)-repmat(kron(ones(1,Nd*Ns),(Nd-1):-1:0),Ne,Ns); J=J';
            S=sparse(J(:),I(:),S(:),Nd*Nc,Nd*Nc);
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
          Nd       = size(co,1);
          Nn       = size(el,2);
          Nq       = size(w, 1);
          Nl       = size(l,2);
          Ns       = size(phi,1);
         
          F    = ofem.matrixarray(zeros(Nd,Ns,Ne));
          elco = reshape(co(:,:,el(:,1:Nl)'),[],Nl,Ne);

            for q=1:Nq
                X = elco*(l(q,:)');
                F = F + f(X)*(w(q)*phi(:,q)');
            end
            F = F*detD;
            I    = repelem(1:Nn,1,Nd); I=I(:);
            I    = Nd*el(:,I) - repmat((Nd-1):-1:0,Ne,Ns); I=I';
            b    = sparse(I(:),1,F(:),Nd*size(co,3),1);
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
            Nd     = size(co ,1);
            Nq     = size(w  ,1);
            Ns     = size(phi,1);
            Nf     = size(faces,1);
            Nn     = size(faces,2);
            
            

           % F      = zeros(Nd*Ns*Nf,1);

            % faceco*l gives global quadrature points => evaluate there
            F      = ofem.matrixarray(zeros(Nd, Ns,Nf));
            faceco = reshape(co(:,:,faces(:,1:Nl)'),[],Nl,Nf);
            
            %faceco  = reshape(permute(reshape(co(faces(:,1:Nl),:),Nf,Nl,Nd),[1,3,2]),Nd*Nf,Nl);

            for q=1:Nq
                X = faceco*(l(q,:)');
                F = F + g(X,normals)*(w(q)*phi(:,q)');
            end

            F    = F*meas;
            I    = repelem(1:Nn,1,Nd); I=I(:);
            I    = Nd*faces(:,I) - repmat((Nd-1):-1:0,Nf,Ns); I=I';
            b    = sparse(I(:),1,F(:),Nd*size(co,3),1);
%             b    = accumarray(faces(:),F);
        end
        
        function b=dirichlet(d,nodes,co)
        %DIRICHLET returns the Dirichlet-originated part of the load vector.
        %
        % b=DIRICHLET(el,co) computes the Dirichlet-originated vector at
        % the specified faces.
        %
%         diri_x    = mesh.dim*dirichlet - 2;
%         diri_y    = mesh.dim*dirichlet - 1;
%         diri_z    = mesh.dim*dirichlet;
            Nc  = size(co ,3);
            dim = size(co,1);
            Dx  =  d(co(1,:,nodes));            
            Dy   = d(co(2,:,nodes));
            Dz   = d(co(3,:,nodes));
            D =squeeze([Dx Dy Dz]);
            D= D';
            pos =[dim*nodes-2 dim*nodes-1 dim*nodes]; 
            b   = sparse(pos(:),1,D(:),dim*Nc,1);
        end
    end

    methods
        %%
        function obj=elastic(mesh,fe, qr)
%             obj@ofem.laplace(mesh,fe);
            obj.mesh  = mesh;
            obj.fe = fe;
            obj.qr = qr;
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
        
            Np  = size(obj.mesh.parts,2);
            Nbd = size(obj.mesh.bd   ,2);
            Nc  = size(obj.mesh.co   ,3);
            Nd  = obj.mesh.dim;
            
            intvol  = 0;
            intface = 0;
            intdiri = 0;
            
            
             %% check input
            if nargin==1
                % elliptic equation with homogenous Neumann boundary
                opt=struct('S',1,'D',0,'M',0,'force',0,'A',1,'b',ones(3,1),'c',1);
            else
                if nargin>3
                    error('ofem:elliptic:TooManyArguments',...
                          'I''m expecting at most one argument!');
                end
                opt=varargin{1};

                if ~isstruct(opt)
                    error('ofem:elliptic:InvalidArgument',...
                          'opt is expected to be a structure.');
                end

                %% stiffness
                if ~isfield(opt,'S')
                    opt.S = 0;
                else
                    % stiffness matrix requested, check for material
                    if isfield(opt, 'A')
                    %if isempty(opt.A)
                        % set to default
                    else
                    end
                    aux.S = cell(Np,1);
                end

                %% damping
                if ~isfield(opt,'D')
                    opt.D=0;
                else
                    aux.D = cell(Np,1);
                end

                %% mass
                if ~isfield(opt,'M')
                    opt.M=0;
                else
                    aux.M = cell(Np,1);
                    if ~isfield(opt,'c')
                        opt.c=1;
                    end
                end

                %% volume force
                if ~isfield(opt,'force')
                    opt.force=[];
                else
                    aux.force = cell(Np,1);
                end

                %% Neumann boundary, pressure
                if ~isfield(opt,'neumann')
                    opt.neumann={};
                else
                    if iscell(opt.neumann)
                        for k=1:numel(opt.neumann)
                            if ~all(isfield(opt.neumann{k},{'data','idx'}))
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.neumann{%d} must contain a ''data'' and a ''idx'' field.',k);
                            end
                            if ~isnumeric(opt.neumann{k}.idx) || ~isscalar(opt.neumann{k}.idx)
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.neumann{%d}.idx must be a scalar.',k);
                            end
                            if opt.neumann{k}.idx>Nbd
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.neumann{%d}.idx exceeds the number of available sidesets.',k);
                            end
                            if isnumeric(opt.neumann{k}.data) && isscalar(opt.neumann{k}.data)
                                val = opt.neumann{k}.data;
                                opt.neumann{k}.data = @(X) val*ones(1,1,size(X,3));
                            elseif isa(opt.neumann{k}.data,'function_handle')
                            else
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.neumann{%d}.data must either be a scalar or a function handle.',k);
                            end
                        end
                    elseif isstruct(opt.neumann) && all(isfield(opt.neumann,{'data','idx'}))
                        neumanndata = opt.neumann.data;
                        neumannidx  = opt.neumann.idx ;
                        
                        if ~isnumeric(neumannidx) || ~isvector(neumannidx)
                            error('ofem:elliptic:InvalidArgument',...
                                  'opt.neumann.idx must be a vector.');
                        end
                        if opt.neumann.idx>Nbd
                            error('ofem:elliptic:InvalidArgument',...
                                  'opt.neumann.idx exceeds the number of available sidesets.');
                        end
                        if isnumeric(neumanndata) && isscalar(neumanndata)
                            val = neumanndata;
                            neumanndata = @(X) val*ones(1,1,size(X,3));
                        elseif isa(neumanndata,'function_handle')
                        else
                            error('ofem:elliptic:InvalidArgument',...
                                  'opt.neumann.data must either be a scalar or a function handle.');
                        end
                        
                        opt.neumann = cell(numel(neumannidx),1);
                        for k=1:numel(neumannidx)
                            opt.neumann{k}.data = neumanndata;
                            opt.neumann{k}.idx  = neumannidx(k);
                        end
                    else
                        error('ofem:elliptic:InvalidArgument',...
                              'opt.neumann must either be a cell array or a structure containing a ''data'' and a ''idx'' field.');
                    end

                    Nneu        = numel(opt.neumann);
                    aux.neumann = cell(Nneu,1);
                end


                %% Dirichlet boundary
                if ~isfield(opt,'dirichlet')
                    opt.dirichlet={};
                else
                    if iscell(opt.dirichlet)
                        for k=1:numel(opt.dirichlet)
                            if ~all(isfield(opt.dirichlet{k},{'data','idx'}))
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.dirichlet{%d} must contain a ''data'' and a ''idx'' field.',k);
                            end
                            if opt.dirichlet{k}.idx>Nbd
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.dirichlet{%d}.idx exceeds the number of available sidesets.',k);
                            end
                            if isnumeric(opt.dirichlet{k}.data) && isscalar(opt.dirichlet{k}.data)
                                val = opt.dirichlet{k}.data;
                                opt.dirichlet{k}.data = @(X) val*ones(1,1,size(X,3));
                            elseif isa(opt.dirichlet{k}.data,'function_handle')
                            else
                                error('ofem:elliptic:InvalidArgument',...
                                      'opt.dirichlet.data must either be a scalar or a function handle.');
                            end
                        end
                    elseif isstruct(opt.dirichlet) && all(isfield(opt.dirichlet,{'data','idx'}))
                        dirichletdata = opt.dirichlet.data;
                        dirichletidx  = opt.dirichlet.idx ;
                        
                        if ~isnumeric(dirichletidx) || ~isvector(dirichletidx)
                            error('ofem:elliptic:InvalidArgument',...
                                  'opt.dirichlet.idx must be a vector.');
                        end
                        if opt.dirichlet.idx>Nbd
                            error('ofem:elliptic:InvalidArgument',...
                                  'opt.dirichlet.idx exceeds the number of available sidesets.');
                        end
                        if isnumeric(dirichletdata) && isscalar(dirichletdata)
                            val = dirichletdata;
                            dirichletdata = @(X) val*ones(1,1,size(X,3));
                        elseif isa(dirichletdata,'function_handle')
                        else
                            error('ofem:elliptic:InvalidArgument',...
                                  'opt.dirichlet.data must either be a scalar or a function handle.');
                        end
                        
                        opt.dirichlet = cell(numel(dirichletidx),1);
                        for k=1:numel(dirichletidx)
                            opt.dirichlet{k}.data = dirichletdata;
                            opt.dirichlet{k}.idx  = dirichletidx{k};
                        end
                    else
                        error('ofem:elliptic:InvalidArgument',...
                              'opt.dirichlet must either be a cell array or a structure containing a ''data'' and a ''idx'' field.');
                    end

                    Ndiri         = numel(opt.dirichlet);
                    aux.dirichlet = cell(Ndiri,1);
                end
            end


            if opt.S||opt.D||opt.M||isempty(opt.force)==0
                intvol=1;
            end

            if isempty(opt.neumann)==0
                intface=1;
            end

            if isempty(opt.dirichlet)==0
                intdiri=1;
            end

            if ~(intvol==1||intface==1||intdiri==1)
                warning('Requested to compute nothing!');
            end
            
% 
%             if nargin==1
%                 opt.S=1;
%                 opt.M=1;
%                 opt.b=0;
%                 opt.N=0;
%             else
%                 opt=varargin{1};
%                 if ~isstruct(opt)
%                     error('opt is expected to be a structure.');
%                 else
% 
%                 if ~isfield(opt,'M'); opt.M=0; end
%                 if ~isfield(opt,'S'); opt.S=0; end
%                 if ~isfield(opt,'b'); opt.b=0; end
%                 if ~isfield(opt,'N'); opt.N=0; end
%                 end
%             end
% 
% %             if opt.M==1
% %                 warning('No mass matrix computation supported so far!');
% %                 opt.M=0; % no stiffness matrix so far!
% %             end
% 
%             if opt.b==1
%                 if ~isfield(opt,'load')
%                     error('Load vector requested, but no function handle given!');
%                 elseif ~isa(opt.load,'function_handle')
%                     error('load is expected to be a function handle!');
%                 end
%             end
%             
%             if opt.N==1
%                 if ~isfield(opt,'neumann')
%                     error('Neumann boundary requested, but no function handle given!');
%                 elseif ~isa(opt.neumann,'function_handle')
%                     error('neumann is expected to be a function handle!');
%                 end
% 
%                 if ~isfield(opt,'neumannidx')
%                     error('You need to specifiy the Neumann boundary index!');
%                 end
%             end

            S  = sparse(Nc*Nd, Nc*Nd);
            M  = sparse(Nc*Nd, Nc*Nd);
            D  = sparse(Nc*Nd, Nc*Nd);
            b  = sparse(Nc*Nd,1);
            
            opt.A=1;
            opt.b=ones(3,1);
           %opt.c=1;

            tic;
%             [DinvT,detD] = obj.mesh.jacobiandata;
%             Np           = size(obj.mesh.parts,2);

%% volume related integration
            if intvol==1
                % volume quad data
                [w,l] = obj.qr.data(0);

                % shape functions related stuff
                pipj  = obj.fe.phiiphij(obj.mesh.dim);
                phi   = obj.fe.phi(l);
                dphi  = obj.fe.dphi(l);
                [DinvT,detD] = obj.mesh.jacobiandata;
                aux.detD   = detD;

                % perform assembly
                for i=1:Np
                    pIdx     = obj.mesh.parts{3,i};
                    elemsLoc = obj.mesh.el(pIdx,:);
                    detDLoc  = detD(:,:,pIdx);
                    DinvTLoc = DinvT(:,:,pIdx);
                    
                    %% handle stiffness matrix
                    if opt.S==1
                        aux.S{i} = obj.stiffness(obj.lambda(i), obj.mu(i),DinvTLoc,detDLoc,dphi,w,elemsLoc,obj.mesh.co);
                        S = S + aux.S{i};
                    end
                    %% handle damping matrix
                    if opt.D==1
                        aux.D{i} = obj.damping(opt.b,DinvTLoc,detDLoc,phi,dphi,w,elemsLoc,obj.mesh.co);
                        D = D + aux.D{i};
                    end
                    %% handle mass matrix
                    if opt.M==1
                        aux.M{i} = obj.mass(opt.c(:,i,pIdx),detDLoc,pipj,elemsLoc,obj.mesh.co);
                        M = M + aux.M{i};
                    end
                end
                %% handle volume force matrix
                if isempty(opt.force)==0
                     aux.force = obj.load(detD,phi,w,l,opt.force{1},obj.mesh.el,obj.mesh.co);
                     b = aux.force;
                end 
            end

            %% surface related integration
            if intface==1
                % surface quad data
                [w,l] = obj.qr.data(1);
                phi   = obj.fe.phi(l);
                aux.neumannidx=zeros(Nneu,1);

                for i=1:Nneu
                    [meas,faces,normals,~] = obj.mesh.neumann(opt.neumann{i}.idx);
                    % [TODO] include the correct function!
                    %neumann(meas,phi,w,l,g,faces,normals,co)
                    aux.neumann{i} = obj.neumann(meas{1},phi,w,l,opt.neumann{i}.data,faces{1},normals{1},obj.mesh.co);
                    aux.neumannidx(i) = opt.neumann{i}.idx;
                    b = b + aux.neumann{i};
                end
            end

            %% Dirichlet data
            if intdiri==1
                aux.dirichletidx=zeros(Ndiri,1);
                for i=1:Ndiri
                    nodes = obj.mesh.dirichlet(opt.dirichlet{i}.idx);
                    aux.dirichlet{i} = obj.dirichlet(opt.dirichlet{i}.data,nodes{1},obj.mesh.co);
                    aux.dirichletidx(i) = opt.dirichlet{i}.idx;
                    b = b - (S+D+M)*aux.dirichlet{i};
                end
            end

            if opt.S
                asm.S=S;
            end
            if opt.D
                asm.D=D;
            end
            if opt.M
                asm.M=M;
            end
            if isempty(opt.force)==0 || isempty(opt.neumann)==0 || isempty(opt.dirichlet)==0
                asm.b=b;
            end


            %% info variable
            info.time2assemble = toc;
        end

        
        function grad=gradu(obj,u)
        %GRADU computes the gradient at DOFs.
        %
        % grad=gradu(u) computes the gradient grad of the FEM solution u at
        % DOFs. grad is a Ndofs by Nd matrix, where Ndofs are the number of
        % nodes and Nd the dimension of the spatial space.
        %
        
        switch obj.fe
            case ofem.finiteelement.P1
                
                switch obj.mesh.dim
                    case 2
                        % d denotes the dimension of polynomial space,
                        % i.e. for P1 elements it is 1
                        d=1;             % degree of finite element space
                        m=(d+2)*(d+3)/2; % polynomial degree of approximant
                        n=(d+3)*(d+4)/2; % degree of point considering for least-squares
                        
                        % statistics and machine learning toolbox
                        co_idx = knnsearch(squeeze(obj.mesh.co)',squeeze(obj.mesh.co)','K',n)';
                        
                        ui = u          (co_idx(:)  );
                        xi = squeeze(obj.mesh.co(1, :, co_idx(:)));
                        yi = squeeze(obj.mesh.co(2, :, co_idx(:)));
                        
                        clear co_idx;
                        
                        Ndof = size(obj.mesh.co,3);
                        
                        A = [ones(Ndof*n,1), xi, yi, xi.*yi, xi.^2, yi.^2];
                        
                        ui_rec = cellfun(@(A,b) A\b      , ...
                            mat2cell(A ,n*ones(1,Ndof),m), ...
                            mat2cell(ui,n*ones(1,Ndof),1), ...
                            'UniformOutput',false);
                        ui_rec = reshape(cell2mat(ui_rec),m,[])';
                        
                        clear A;
                        x = squeeze(obj.mesh.co(1, :,:));
                        y = squeeze(obj.mesh.co(2,:, :));
                        
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
                       co_idx = knnsearch(squeeze(obj.mesh.co)',squeeze(obj.mesh.co)','K',n)';
                        
                        %                 ui = u (co_idx(:)  );
                        %                 xi = co(co_idx(:),1);
                        %                 yi = co(co_idx(:),2);
                        %                 zi = co(co_idx(:),3);
                        
                        
                        
                        
                        ui = u          (co_idx(:)  );
                        xi = squeeze(obj.mesh.co(1, :, co_idx(:)));
                        yi = squeeze(obj.mesh.co(2, :, co_idx(:)));
                        zi = squeeze(obj.mesh.co(3, :, co_idx(:)));
                        
                        
                        clear co_idx;
                        
                        Ndof = size(obj.mesh.co,3);
                        
                        A = [ones(Ndof*n,1), xi, yi, zi, xi.*yi ,xi.*zi, yi.*zi, xi.^2, yi.^2, zi.^2];
                        
                        ui_rec = cellfun(@(A,b) A\b      , ...
                            mat2cell(A ,n*ones(1,Ndof),m), ...
                            mat2cell(ui,n*ones(1,Ndof),1), ...
                            'UniformOutput',false);
                        ui_rec = reshape(cell2mat(ui_rec),m,[])';
                        
                        clear A;
                        
                        x = squeeze(obj.mesh.co(1, :,:));
                        y = squeeze(obj.mesh.co(2,:, :));
                        z = squeeze(obj.mesh.co(3, :,:));
                        
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
            
            F = zeros(obj.mesh.dim*size(obj.mesh.co,3),obj.mesh.dim);
            for i=1:obj.mesh.dim
                F(i:obj.mesh.dim:end,:) = obj.gradu(u(:,i));
            end
            for i=1:obj.mesh.dim
                F(i:obj.mesh.dim:end,i) = F(i:obj.mesh.dim:end,i)+1;
            end
        end
        
        function [E,e,s,S] = StrainStress(obj,u)
            
            co     = obj.mesh.co;
            dim    = obj.mesh.dim;
            F      = obj.defGrad(u);
            Np     = size(obj.mesh.parts,2);
            s      = cell(Np,1);
            e      = cell(Np,1);
            S      = zeros(size(u,1), 6);
            
            %% Green's Stain Tensor computation ( E =1/2*(F^T F - I))
            K= permute(reshape(F,dim,size(co,3),dim),[1 3 2]);
            I=repmat(1:dim,dim,1);
            J=I';
            EE=reshape(dot(K(I,:,:), K(J,:,:),2),dim,dim,[]);
            %this is the strain now, but in the wrong format so we have to change this
            E=0.5*(reshape(permute(EE,[1,3,2]),[],dim)-repmat(eye(dim,dim),...
                size(co,3),1));
            K= permute(reshape(E,dim,size(co,3),dim),[1 3 2]);
            

            E = reshape(K,dim*dim, size(u,1))';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %    Quick and Dirty Workaround for Stress tensor computation
            %    For more than one material, only the first material
            %    appearing in the list is considered for the stress
            %    computation at the interface
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j=1:Np
                pIdx = obj.mesh.parts{3,j};
                elems{j,:} = unique(obj.mesh.el(pIdx,:));
            end
            [interElems, ~] = intersectAll(elems);
           
            switch dim
                case 3
                     %isotropic stress tensor: only the idependen components need to be
                     %computed, this happens by using only the column components
                     %[1 5 9 6 3 2] (i.e. e_{11} e_{22} e_{33} e_{23} e_{13} e_{12}) of every row
                    for i =1:Np
                        if(i>1)
                            for rowI=1:(i-1)
                                nIdx =setdiff(elems{i,:}, interElems{rowI, i});
                            end
                        else
                            nIdx = elems{i,:};
                        end
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
                        S(nIdx,:) =S(nIdx,:) +s{i}'; 
                    end
                case 2
                    warning('MATLAB:ofem:elastic:StrainStress:WordsOfWisdom',...
                            'You are calculating stresses in 2D, only plane stress is implemented so far, I will compute that')
                    %isotropic stress tensor: only the idependen components need to be
                    %computed, this happens by using only the column components
                    %[1 4 2] (i.e. e_{11} e_{22} e_{12}) of every row
                    for i =1:Np   
                        if(i>1)
                            for rowI=1:(i-1)
                                nIdx =setdiff(elems{i,:}, interElems{rowI, i});
                            end
                        else
                            nIdx = elems{i,:};
                        end
%                         pIdx     = obj.mesh.parts{3,i};
%                         elemsLoc = obj.mesh.el(pIdx,:);
%                         nIdx     = unique(elemsLoc);
                        uPart    = u(nIdx, :);
                        e{i} = repmat([ 1, 1,2],size(uPart,1), 1).*...
                            E(nIdx, [1,4,2]);
                        C = [2*obj.mu(i) + obj.lambda(i),  obj.lambda(i), 0;
                            obj.lambda(i),2*obj.mu(i)+obj.lambda(i),  0;
                             0,            0,                 obj.mu(i)];
                        
                        s{i}=C*e{i}';
                        S(nIdx,:) =S(nIdx,:) +s{i}'; 
                    end
                case 1
            end
            
        end
    end
end