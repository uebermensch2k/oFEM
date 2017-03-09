classdef gaussianquadrature
    properties(Access=protected)
        dim;
        fe;
    end

    methods(Access=protected,Static)
        %%
        function [w,l] = quaddata(fe,dim)

            switch fe
                case ofem.finiteelement.P1
                    w=1/factorial(dim);
                    l=repmat(1/(dim+1),1,dim+1);

                case ofem.finiteelement.P2
                    switch dim
                        case 0
                            w = 1;
                            l = 1;
                            
                        case 1
                            w = [1; 1]/2;
                            l = 0.5+[-1; 1]/(2*sqrt(3));
                            l = [1-l, l];

                        case 2
                            w = [1; 1; 1]/6;
                            l = [1/6 1/6; 2/3 1/6; 1/6 2/3];
                            l = [1-sum(l,2),l];

                        case 3
                            w = [-2/15; 3/40; 3/40; 3/40; 3/40];
                            l = [ 1/4 1/4 1/4; ...
                                  1/6 1/6 1/6; ...
                                  1/6 1/6 1/2; ...
                                  1/6 1/2 1/6; ...
                                  1/2 1/6 1/6];
                            l = [1-sum(l,2),l];
                    end
                    
                case ofem.finiteelement.Q1
                    error('ofem:gaussianquadrature:quaddata:NotImplemented',...
                          'Quassian quadrature for Q1 elements not implemented yet!');
                    
                case ofem.finiteelement.Q2
                    error('ofem:gaussianquadrature:quaddata:NotImplemented',...
                          'Quassian quadrature for Q2 elements not implemented yet!');
                    
                case ofem.finiteelement.NE0P
                    error('ofem:gaussianquadrature:quaddata:NotImplemented',...
                          'Quassian quadrature for NE0P elements not implemented yet!');

                case ofem.finiteelement.NE1P
                    error('ofem:gaussianquadrature:quaddata:NotImplemented',...
                          'Quassian quadrature for NE1P elements not implemented yet!');

                case ofem.finiteelement.NE2P
                    error('ofem:gaussianquadrature:quaddata:NotImplemented',...
                          'Quassian quadrature for NE2P elements not implemented yet!');
            end
        end
    end

    methods
        %%
        function obj=gaussianquadrature(mesh,fe)
        %gaussianquadrature construct gaussianquadrature
        %
        % see also ofem.mesh, ofem.finiteelement
        %
            obj.dim = mesh.dim;
            obj.fe  = fe;
        end

        %%
        function [w,l] = data(obj,codim)
        %data returns quadrature points and weights.
        %
        % [w,l]=data(fe,codim) returns per row the barycentric coordinates
        % l of the quadrature points and weights w sufficient for the
        % specified element and codimension. The quadrature rule is a
        % Gaussian quadrature. The quadarture points l yield barycentric
        % coordinates as they are used by ofem.finiteelement.phi and
        % ofem.finiteelement.dphi.
        %
        % see also ofem.finiteelement.phi, ofem.finiteelement.dphi
        %
            edim = obj.dim-codim;
            if edim<0
                error('ofem:gaussianquadrature:data',...
                      'Codimension cannot exceed dimension of embedded space!');
            end
            [w,l]=quaddata(obj.fe,obj.dim-codim);
        end
    end
end

