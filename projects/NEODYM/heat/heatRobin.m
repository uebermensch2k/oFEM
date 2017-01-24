function T = heatRobin(xpts,tpts,f,T0,varargin)
%heatRobin computes solution of heat equation
%
% T=
%
% U=heatRobin(x,t,f,T0,coeffs), given an array of function handles f, computes
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
% By default, a unit square is assumed as the domain and oak as medium with
% following characteristics:
%   - c_P = 0.17 J/(kg K)  (specific heat)
%   - rho = 750 kg/m^3     (density)
%   - k   = 2400 W/(m K) (thermal conductivity)
% and ambient temperature T_a=298.15 K (i.e., 25 degree C).
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
    addRequired(p,'t',@(x) validateattributes(x,{'numeric'},{'nonempty','real','vector','increasing','>=',0}));
    addRequired(p,'f',@(x) validateattributes(x,{'function_handle'},{'nonempty','vector'}));
    addRequired(p,'coeffs',@(x) validateattributes(x,{'struct'},{'nonempty'}));
    addOptional(p,'domain',struct('LL',[0,0],'UR',[1,1],'h_max',1/20),@(x) validateattributes(x,{'struct'},{'nonempty'}));
    parse(p,xpts,tpts,f,varargin{:});

    domain=p.Results.domain;
    coeffs=p.Results.coeffs;

    if ~all(isfield(coeffs,{'v','h'}))
        error('heatRobin:argchk','Input argument ''coeffs'' need to be a structure with containing at least the fields ''v'', ''h''');
    else
        if ~isfield(coeffs,'c_P')
            coeffs.c_P = 0.17;
        end
        if ~isfield(coeffs,'rho')
            coeffs.rho = 750;
        end
        if ~isfield(coeffs,'k')
            coeffs.k = 2400;
        end
        if ~isfield(coeffs,'T_a')
            coeffs.T_a = 298.15;
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

%     info = mesh.info;
%     dx = info.edge.max;

    % assign finit element
    fe = ofem.finiteelement.P1;
    mesh.assign_fe(fe);

    % differential operator
    op = ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));

    % assemble

    % Be careful, the mass matrix comes from the time term. Thus do not
    % include mass in Dirichlet conditions.
    opt = struct('M',1,'c',1);
    [asm,~,~] = op.assemble(opt);
    M       = asm.M;
    clear opt;

    
    % stiffness matrix
    opt.S     = 1;
    opt.A     = coeffs.k/(coeffs.c_P*coeffs.rho);
    % damping matrix
    opt.D     = 1; % compute damping
    opt.b     = coeffs.v;
    % Robin conditions
    opt.robin.idx   = 1;
    opt.robin.alpha = coeffs.h;
    opt.robin.f     = coeffs.h*coeffs.T_a;
    % force
    opt.force       = @(x) f(x)./(coeffs.c_P*coeffs.rho);
    % do the job
    [asm,~,~] = op.assemble(opt);

    S       = asm.S;
    D       = asm.D;
    M_robin = asm.M_robin;
    b       = asm.b;
    DOFs    = asm.DOFs; % no Dirichlet boundary, therefore no DOFs needed
    A       = S(DOFs,DOFs)+D(DOFs,DOFs)+M_robin(DOFs,DOFs);
%     clear S D M_robin;

    % initial condition
    T0 = squeeze(permute(double(T0(mesh.co)),[3,1,2]));

    % solve
%     if max(norm(coeffs.v))~=0
%         odeopt = odeset('Mass',M(DOFs,DOFs),'Stats','on','NonNegative',1:size(M,1),'MaxStep',dx^2/(max(norm(coeffs.v))));
%     else
%         odeopt = odeset('Mass',M(DOFs,DOFs),'Stats','on','NonNegative',1:size(M,1));
%     end
    odeopt = odeset('Mass',M(DOFs,DOFs),'Stats','on','NonNegative',1:size(M,1));
    sol    = ode45(@(t,x) b(DOFs)-A*x,tpts,T0(DOFs),odeopt);
    clear M b;
    T_t    = deval(sol,tpts);
    clear sol DOFs;

    TR     = triangulation(mesh.el,permute(double(mesh.co),[3,1,2]));
    [i,l]  = pointLocation(TR,xpts);
    clear TR;
    T      = reshape(T_t(mesh.el(i,:),:),[],size(mesh.el,2),numel(tpts));
    T      = reshape(dot(T,repmat(fe.phi(l)',1,1,numel(tpts)),2),[],numel(tpts));

    visualize=1;
    if visualize
        figure
        view(50,45);
        xlabel('x-axis [m]');
        ylabel('y-axis [m]');
        zlabel('T [K]');
        hold on;
        for i=1:numel(tpts)
            cla;
            zlim([min(T_t(:,i))-eps(min(T_t(:,i))),max(T_t(:,i))+eps(max(T_t(:,i)))]);
            trimesh(mesh.el,double(mesh.co(1,1,:)),double(mesh.co(2,1,:)),T_t(:,i));
            plot3(xpts(:,1),xpts(:,2),T(:,i),'r*');
            title(sprintf('t=%f',tpts(i)));
            pause(0.1);
        end
    end
    
    
end