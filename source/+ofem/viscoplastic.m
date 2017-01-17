classdef viscoplastic < ofem.elastic
    properties(Access=protected)
%         mesh;
%         felem;
%         lambda; %Lame's 1st constant
%         mu;     %Lame's 2nd constant
        nu;     %plastic viscosity
        sigY;   %Yield Stress

        theta;
        dt;
    end
    
    methods(Access=protected,Static)
        
        function [M,Q] = stiffness(DinvT,detD,dphi,el,co,...
                u1,u0,sig0,sigY,lam,mu,nu,theta,dt,w)
        % NEWTON_MATRICES returns the matrices M and F used for the
        % Newton iteration in the case of perfect viscoplasticity
        % without hardening.
        % F denotes the function resulting from the weak formulation of
        % the viscoplastic problem for which we seek the root. M is the
        % Jacobian of F with respect to the solution vector.
        % For a detailed discussion of the used model see:
        %   Elastoviscoplastic Finite Element analysis in 100 lines of
        %   Matlab
        %   C. Carstensen and R. Klose
        %   J. Numer. Math., Vol. 10, No.3 pp.157-192 (2002)
        %
            Ns = size(dphi,1);
            Nd = size(dphi,2);
            Nq = size(dphi,3);
            Ne = size(el,1);
            Nc = size(co,1);

            C1 = lam+2*mu/Nd;
            C2 = nu/(nu/(2*mu)+theta*dt);
            C3 = theta*dt*sigY/(nu/(2*mu)+theta*dt);

            e0=1/(9*lam+6*mu)*ofem.tr(sig0,Nd)*reshape(eye(Nd,Nd),1,[])+1/(2*mu)*ofem.dev(sig0,Nd);

            IF = 3*repelem(el,1,Nd)-kron(ones(Ne,1),repmat((Nd-1):-1:0,1,Ns));

            udiff = reshape(u1(IF')-u0(IF'), Nd*Ns, 1, Ne);

            M=zeros(Nd*Ns,Nd*Ns,Ne);
            Q=zeros(Nd*Ns,1    ,Ne);
            for k=1:Nq
                % DinvT is Nd*NexNd matrix
                % dphi is NsxNdxNl matrix
                % dphi(:,:,k)*DinvT' is NsxNd*Ne
                globdphi=reshape(dphi(:,:,k)*DinvT',Ns,Nd,Ne);

                Deta1 = zeros(Ns*Nd,Nd*Nd,Ne);
                Deta2 = zeros(Ns*Nd,Nd*Nd,Ne);

                for i=1:Nd
                    Deta1(i:Nd:end,Nd*(i-1)+1:i*Nd,:) = globdphi;
                    Deta2(i:Nd:end,i:Nd:end       ,:) = globdphi;
                end

                eps=(Deta1+Deta2)/2;

                v=ofem.mult(permute(udiff,[2,1,3]),eps)+permute(e0,[3,2,1]);

                % elastic
                C5=2*mu*ones(1,1,Ne);
                C6=zeros(Ns*Nd,1,Ne);

                % plastic
                devv=ofem.dev(v,Nd);
                ndevv=sqrt(dot(devv,devv,2));
                plastic_idx=squeeze(ndevv>sigY/(2*mu));

                if any(plastic_idx)
                    fprintf('Plastic phase for %d indices\n',sum(plastic_idx));

                    devv  = devv    (:,:,plastic_idx);
                    ndevv = ndevv   (:,:,plastic_idx);

                    deve  = ofem.dev(eps(:,:,plastic_idx),Nd);

                    % per element: C5=C2+C3/norm(dev3(v));
                    C5(:,:,plastic_idx) = C2+C3./ndevv;

                    % per element: C6=C3/norm(dev3(v))^3*(dev3(eps)*dev3(v)');
                    C6(:,:,plastic_idx) = ofem.mult(C3./ndevv.^3,ofem.mult(deve,permute(devv,[2,1,3])));
                end

                % per element: M=M0+T*(C1*tr3(eps)*tr3(eps)'+C5*dev3(eps)*eps'-C6*dev3(v)*eps');
                M=M+w(k)*(...
                           C1*ofem.mult(ofem.tr(eps,Nd),permute(ofem.tr(eps,Nd),[2,1,3])                  ) ...
                          +   ofem.mult(C5             ,ofem.mult(ofem.dev(eps,Nd), permute(eps,[2,1,3])) ) ...
                          -   ofem.mult(C6             ,ofem.mult(ofem.dev(v  ,Nd), permute(eps,[2,1,3])) ) ...
                          );

                % per element: F=F0+T*(C1*tr3(v)*tr3(eps)+C5*eps*dev3(v)');
                Q=Q+w(k)*(...
                           C1*ofem.mult(ofem.tr(v,Nd), ofem.tr(eps,Nd)                                ) ...
                          +   ofem.mult(C5           , ofem.mult(eps,permute(ofem.dev(v,Nd),[2,1,3])) ) ...
                         );

            end
            
            M=ofem.mult(M,permute(detD,[2,3,1]));
            Q=ofem.mult(Q,permute(detD,[2,3,1]));
            
            I  = repmat (IF,1,Nd*Ns)';
            J  = repelem(IF,1,Nd*Ns)';
            IF = IF';
            
            M = sparse(I(:),J(:),M(:),Nd*Nc,Nd*Nc);
            Q = sparse(IF(:),1,Q(:),Nd*Nc,1);
        end

        %%
        function [sigma] = tension(DinvT,dphi,el,co,...
                u1,u0,sig0,sigY,lam,mu,nu,theta,dt,w)
            % TENSION returns the stress tensor computed by disassambling
            %
            % Elastoviscoplastic Finite Element analysis in 100 lines of
            % Matlab
            % C. Carstensen and R. Klose
            % J. Numer. Math., Vol. 10, No.3 pp.157-192 (2002)
            
            
            Ns = size(dphi,1);
            Nd = size(dphi,2);
            Nq = size(dphi,3);
            Ne = size(el,1);
%             Nc = size(co,1);

            C1 = lam+2*mu/Nd;
            C2 = nu/(nu/(2*mu)+theta*dt);
            C3 = theta*dt*sigY/(nu/(2*mu)+theta*dt);

            e0=1/(9*lam+6*mu)*ofem.tr(sig0,Nd)*reshape(eye(Nd,Nd),1,[])+1/(2*mu)*ofem.dev(sig0,Nd);
            
            IF = 3*repelem(el,1,Nd)-kron(ones(Ne,1),repmat((Nd-1):-1:0,1,Ns));
            
            udiff = reshape (u1(IF')-u0(IF'), Nd*Ns, 1, Ne);
            sigma=zeros(1,Nd*Nd,Ne);
            for k=1:Nq
                % DinvT is Nd*NexNd matrix
                % dphi is NsxNdxNl matrix
                % dphi(:,:,k)*DinvT' is NsxNd*Ne
                globdphi=reshape(dphi(:,:,k)*DinvT',Ns,Nd,Ne);
                
                Deta1 = zeros(Ns*Nd,Nd*Nd,Ne);
                Deta2 = zeros(Ns*Nd,Nd*Nd,Ne);
                
                for i=1:Nd
                    Deta1(i:Nd:end,Nd*(i-1)+1:i*Nd,:) = globdphi;
                    Deta2(i:Nd:end,i:Nd:end       ,:) = globdphi;
                end
                
                eps=(Deta1+Deta2)/2;
                
                v=ofem.mult(permute(udiff,[2,1,3]),eps)+permute(e0,[3,2,1]);

                dev3v  = ofem.dev(v,Nd);
                ndev3v = sqrt(dot(dev3v,dev3v,2));

                % elastic
                C5=2*mu*ones(1,1,Ne);

                % plastic
                plastic_idx=squeeze(ndev3v>sigY/(2*mu));

                if any(plastic_idx)
%                     dev3v  = dev3v (:,:,plastic_idx);
%                     ndev3v = ndev3v(:,:,plastic_idx);

                    C5(:,:,plastic_idx) = C2+C3./ndev3v(:,:,plastic_idx);
                end
                
                
%                 sigma = sigma + w(k)*(ofem.mult(C1*ofem.tr(v,Nd),reshape(eye(Nd,Nd),1,[]))+ofem.mult(C5,dev3v));
                sigma = sigma + (ofem.mult(C1*ofem.tr(v,Nd),reshape(eye(Nd,Nd),1,[]))+ofem.mult(C5,dev3v));
            end
            sigma=sigma/Nq;

            sigma=squeeze(sigma)';
        end
    end

    methods
        %%
        function obj=viscoplastic(mesh,felem)
        %VISCOLASTIC construct the object from a OFEM.MESH and a
        %OFEM.FINITEELEMENT class.
        %
        % see also OFEM.MESH, OFEM.FINITEELEMENT
        %
            obj@ofem.elastic(mesh,felem);
        end

        %%
        function setmaterial(obj,lambda,mu,nu,sigY)
            setmaterial@ofem.elastic(obj,lambda,mu);
            obj.nu     = nu    ; % plastic viscosity
            obj.sigY   = sigY  ; % Yield Stress
        end

        %%
        function settime(obj,theta,dt)
            obj.theta = theta;
            obj.dt    = dt   ;
        end

        %%
        function [asm,info,aux] = assemble(obj,u1,u0,sig0,varargin)
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
        % Note that the DOFs are Nd x Nc. The solution is therefore
        % orderer as [u1T1;u2T1;u3T1;u1T2;u2T2;u3T2;...]. 
        %

            asm=[];

            if nargin<4
                error('MATLAB:ofem:viscoplastic:assemble:ArgChck',...
                    'Number of arguments is at least four');
            elseif nargin==4
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

            % we have no mass matrix at all
            opt.M=0;

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
%             if opt.M==1
%                 M    = cell(Np,1);
%                 pipj = obj.felem.phiiphij(obj.mesh.dim);
%             end

            if opt.S==1 || opt.b==1
               [w,l] = obj.felem.quaddata(obj.mesh.dim);
            end

            %% stiffness
            if opt.S==1
                S    = cell(Np,1);
                Q    = cell(Np,1);
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
            if opt.S==1 %|| opt.M==1
                for i =1:Np
                    pIdx     = obj.mesh.parts{3,i};
                    elemsLoc = obj.mesh.el(pIdx,:);
                    detDLoc  = detD(pIdx);
                    
                    if opt.S==1
                        Didx        = false(1,length(pIdx)); Didx(pIdx)=true;
                        Didx        = repmat(Didx,obj.mesh.dim,1);
                        DinvTLoc    = DinvT(Didx,:);
                        [S{i},Q{i}] = obj.stiffness(DinvTLoc,detDLoc,dphi,elemsLoc,obj.mesh.co,u1,u0,sig0,obj.sigY,obj.lambda,obj.mu,obj.nu,obj.theta,obj.dt,w);
                    end
%                     if opt.M==1
%                         M{i}     = obj.mass(pipj,detDLoc,elemsLoc,obj.mesh.co);
%                     end
                end
            end
            %% auxiliary variable
            info.time2assemble = toc;

            %% computed structure
            if opt.S==1; asm.S=S; asm.Q=Q; end
%             if opt.M==1; asm.M=M; end
            if opt.b==1 || opt.N==1; asm.b=b; end

            aux.detD=detD;
        end

        function sigma = stress(obj,u1,u0,sig0)
            tic;
%             [DinvT,detD] = obj.mesh.jacobiandata;
            [DinvT,~] = obj.mesh.jacobiandata;
            Np           = size(obj.mesh.parts,2);

            [w,l] = obj.felem.quaddata(obj.mesh.dim);
            dphi  = obj.felem.dphi(l);

            sigma    = cell(Np,1);
            
            for i =1:Np
                pIdx     = obj.mesh.parts{3,i};
                elemsLoc = obj.mesh.el(pIdx,:);
%                 detDLoc  = detD(pIdx);
                
                Didx        = false(1,length(pIdx)); Didx(pIdx)=true;
                Didx        = repmat(Didx,obj.mesh.dim,1);
                DinvTLoc    = DinvT(Didx,:);
                
                sigma{i}    = obj.tension(DinvTLoc,dphi,elemsLoc,obj.mesh.co,u1,u0,sig0,obj.sigY,obj.lambda,obj.mu,obj.nu,obj.theta,obj.dt,w);
            end
            %% auxiliary variable
            info.time2assemble = toc;
        end
    end
end