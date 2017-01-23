function T_xt = heatRobin(xpts,tpts,f,T0,varargin)
%heatRobin computes solution of heat equation
%
% U=heatRobin(x,t,f,coeffs), given an array of function handles f, computes
% for each function handle f_i the solution U at spatial points x and time
% instances t. The underlying PDE is a convection-diffusion equation
%
%   c_P rho \frac{\partial T}{\partial t}-\nabla\cdot(k\nabla T)+v\cdot\nabla T=s
%
% with boundary conditions
%
%   \frac{\partial T}{\partia n} = h(T-T_a),
%
% with specific heat c_P, dnesity rho, thermal conductivity k, wind
% velocity v, heat transfer coefficient h and ambient temperature T_a.
%
% By default, a unit square is assumed as the domain and air as medium with
% following characteristics:
%   - c_P = 1005 J/(kg K)  (specific heat)
%   - rho = 1.2 kg/m^3     (density)
%   - k   = 0.0257 W/(m K) (thermal conductivity)
% and ambient temperature T_a=293.15 K (i.e., 20 degree C).
% However, the wind velocity v and heat transfer coefficient alpha need to
% explicitly be specified. All coeffients are assumed as fields in the
% structure coeffs.
%
% U=heatRobin(x,t,f,domain), additionaly changes the domain to a rectangle
% with lower left corner domain.LL and upper right corner domain.UR. The
% length of the longest edge can be specified by h_max and defaults to
% diagonal length/20.
%

    p=inputParser;
    addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'nonempty','real','2d'}));
    addRequired(p,'t',@(x) validateattributes(x,{'numeric'},{'nonempty','real','vector','increasing'}));
    addRequired(p,'f',@(x) validateattributes(x,{'function_handle'},{'nonempty','vector'}));
    addRequired(p,'coeffs',@(x) validateattributes(x,{'struct'},{'nonempty'}));
    addOptional(p,'domain',struct('LL',[0,0],'UR',[1,1],'h_max',1/20),@(x) validateattributes(x,{'struct'},{'nonempty'}));
%     addOptional(p,'coeffs',struct('c_P',1005,'rho',1.2,'D',0.0257,'a',[5.1;0],'h',1,'T_a',293.15),@(x) validateattributes(x,{'struct'},{'nonempty'}));
    parse(p,xpts,tpts,f,varargin{:});

    domain=p.Results.domain;
    coeffs=p.Results.coeffs;

    if ~all(isfield(coeffs,{'v','h'}))
        error('heatRobin:argchk','Input argument ''coeffs'' need to be a structure with containing at least the fields ''v'', ''h''');
    else
        if ~isfield(coeffs,'c_P')
            coeffs.c_P = 1005;
        end
        if ~isfield(coeffs,'rho')
            coeffs.rho = 1.2;
        end
        if ~isfield(coeffs,'k')
            coeffs.k = 0.0257;
        end
        if ~isfield(coeffs,'T_a')
            coeffs.T_a = 293.15;
        end
    end

    if ~all(isfield(domain,{'LL','UR'}))
        error('heatRobin:argchk','Input argument domain need to be a structure with fields ''LL'' and ''UR''');
    else
        if numel(domain.LL)~=size(xpts,2)
            error('heatRobin:argchk','Computational domain is %d-d, but spatial query points are %d-d.',numel(domain.LL),size(xpts,2));
        end
    end

    if ~isfield(domain,'h_max')
        domain.h_max = norm(domain.UR-domain.LL)/20;
    end

    if numel(tpts)==1
        tpts=[0,tpts];
    end

    % check if points lie inside domain
    y=xpts-repmat(domain.LL,size(xpts,1),1);
    UR=domain.UR-domain.LL;
    if ~all(y>=0 & y<=repmat(UR,size(xpts,1),1))
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

    % Be careful, the mass matrix comes from the time term. Thus do not
    % include mass in Dirichlet conditions.
    opt = struct('M',1,'c',coeffs.c_P*coeffs.rho);
    [asm,~,~] = op.assemble(opt);
    M       = asm.M;
    clear opt;

    
    % stiffness matrix
    opt.S     = 1;
    opt.A     = coeffs.k;
    % damping matrix
    opt.D     = 1; % compute damping
    opt.b     = coeffs.v;
    % Robin conditions
    opt.robin.idx   = 1;
    opt.robin.alpha = coeffs.h;
    opt.robin.f     = coeffs.h*coeffs.T_a;
    % force
    opt.force       = f;
    % do the job
    [asm,~,~] = op.assemble(opt);

    S       = asm.S;
    D       = asm.D;
    M_robin = asm.M_robin;
    b       = asm.b;
    DOFs    = asm.DOFs; % no Dirichlet boundary, therefore no DOFs needed

    % initial condition
    T0 = squeeze(permute(double(T0(mesh.co)),[3,1,2]));

    odeopt = odeset('Mass',M,'Stats','on','NonNegative',1:size(M,1));
    A      = S(DOFs,DOFs)+D(DOFs,DOFs)+M_robin(DOFs,DOFs);
    sol    = ode45(@(t,x) b(DOFs)-A*x,tpts,T0(DOFs),odeopt);
    T_xt   = deval(sol,tpts);

    visualize=1;
    if visualize
        figure
        for i=1:numel(tpts)
            trimesh(mesh.el,double(mesh.co(1,1,:)),double(mesh.co(2,1,:)),T_xt(:,i));
            title(sprintf('t=%f',tpts(i)));
            pause(0.1);
        end
    end
    
    TR=triangulation(mesh.el,permute(double(mesh.co),[3,1,2]));
    [elidx,lambda] = pointLocation(TR,xpts);
    T_xt = reshape(T_xt(mesh.el(elidx,:),:),[],size(mesh.el,2),numel(tpts));
    T_xt = squeeze(dot(T_xt,repmat(fe.phi(lambda)',1,1,numel(tpts)),2));    
end