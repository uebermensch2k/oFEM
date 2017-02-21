function [u,du,ddu] = initialconditions(mesh)

Nc = size(mesh.co,1);

u   = zeros(Nc,1);
du  = zeros(Nc,1);
ddu = zeros(Nc,1);

radius = 5e-2;

idx = dot(mesh.co,mesh.co,2)<radius^2;
if ~isempty(idx)
    du(idx) = 1e-2; % velocity = 10cm/s
%     u(idx) = -1e-2; % displacement = -10cm
end

end