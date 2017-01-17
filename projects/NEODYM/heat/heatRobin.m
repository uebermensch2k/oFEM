function T = heatRobin(x,t,f,varargin)
%heatRobin computes solution of heat equation
%
% U=heatRobin(x,t,f), given an array of function handles f, computes for
% each function handle f_i the solution U at spatial points x and time
% instances t. The underlying PDE is a convection-diffusion equation
%
%   \frac{\partial T}{\partial t}-\nabla\cdot(D\nabla T)+a\cdot\nabla T=s
%
% with boundary conditions
%
%   \frac{\partial T}{\partia n} = h(T-T_a)
%
% with ambient temperature T_a. By default, the diffusion coefficient D is
% set to 1 m^2/s, the wind velocity a is set to [1,0] m/s, and the ambient
% temperture T_a is set to 281.15 K. The domain is the unit square [0,1]^2.
%
% U=heatRobin(x,t,f,domain), additionaly changes the domain to a rectangle
% with lower left corner domain.LL and upper right corner domain.UR.
%
% U=heatRobin(x,t,f,domain,coeffs), additionaly changes the coefficients to
% coeffs.D, coeffs.v and coeffs.T_a, respectively.

    p=inputParser;
    addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'nonempty','real','2d'}));
    addRequired(p,'t',@(x) validateattributes(x,{'numeric'},{'nonempty','real','1d','increasing'}));
    addRequired(p,'f',@(x) validateattributes(x,{'function_handle'},{'nonempty','1d'}));
    addOptional(p,'domain',struct('LL',[0,0],'UR',[1,1],'h_max',1/100),@(x) validateattributes(x,{'struct'},{'nonempty'}));
    addOptional(p,'coeffs',struct('c_P',1,'rho',1,'D',1,'a',[1;0],'h',1,'T_a',281.15),@(x) validateattributes(x,{'struct'},{'nonempty'}));
    parse(p,x,t,f,varargin{:});

    domain=p.Results.domain;
    coeffs=p.Results.coeffs;

    if ~all(isfield(domain,{'LL','UR'}))
        error('heatRobin:argchk','Input argument domain need to be a structure with fields ''LL'' and ''UR''');
    else
        if numel(domain.LL)~=size(x,2)
            error('heatRobin:argchk','Computational domain is %d-d, but spatial query points are %d-d.',numel(domain.LL),size(x,2));
        end
    end

    if ~all(isfield(coeffs,{'c_P','rho','D','a','h','T_a'}))
        error('heatRobin:argchk','Input argument domain need to be a structure with fields ''c_P'', ''rho'', ''D'', ''a'', ''h'' and ''T_a''');
    end

    if ~isfield(domain,'h_max')
        domain.h_max = norm(domain.UR-domain.LL)/100;
    end

    % check if points lie inside domain
    y=x-repmat(domain.LL,size(x,1),1);
    UR=domain.UR-domain.LL;
    if ~all(y>=0 & y<=repmat(UR,size(x,1),1))
%         if any(strcmp(p.UsingDefaults,'domain'))
%             % resize domain
%             domain.LL=[min(x,[],1),min(x,[],2)];
%             domain.UR=[max(x,[],1),max(x,[],2)];
%             offset=0.1*(domain.UR-domain.LL);
%             domain.LL=domain.LL-offset*[1,1];
%             domain.UR=domain.UR+offset*[1,1];
%         else
%             error('heatRobin:argchk','Given points must lie inside ');
%         end
        error('heatRobin:argchk','Given points must lie inside ');
    end
    clear y UR;

    % create mesh
    mesh=ofem.mesh;
    mesh.hypercube(domain.LL,domain.UR-domain.LL);
    k=ceil(log2(norm(domain.UR-domain.LL)/domain.h_max));
    for i=1:k
        mesh.uniform_refine();
    end

    % assign finit element
    fe = ofem.finiteelement.P1;
    mesh.assign_fe(fe);

    % differential operator
    op = ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));

    % assemble
    opt.S     = 1; % compute stiffness
    opt.D     = 1; % compute damping
    opt.M     = 1; % mass matrix
    opt.A     = coeffs.D;
    opt.b     = coeffs.a;
    opt.c     = coeffs.c_P*coeffs.rho;
    opt.robin = struct('idx',1,'alpha',-coeffs.h,'f',-coeffs.h*coeffs.T_a);

    [asm,~,~] = op.assemble(opt);
    asm
    info
    aux

%     ode45(@(x)
end