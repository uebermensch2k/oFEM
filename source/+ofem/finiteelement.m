

classdef finiteelement
    enumeration
        P1, P2, Q1, Q2, NE0P, NE1P, NE2P
    end

    methods
        %%
        function [w,l] = quaddata(obj,dim)
        %QUADDATA returns quadrature points and weights.
        %
        % [w,l]=QUADDATA(dim) returns per row the barycentric coordinates l
        % of the quadrature points and weights w sufficient for the
        % specified element and space dimension given through dim. The
        % quadrature rule is a Gaussian quadrature. The quadarture points l
        % yield barycentric coordinates as they are used by
        % phi and dphi.
        %
        % see also OFEM.FINITEELEMENT.PHI, OFEM.FINITEELEMENT.DPHI
        %

            switch obj
                case ofem.finiteelement.P1
                    w=1/factorial(dim);
                    l=repmat(1/(dim+1),1,dim+1);

                case ofem.finiteelement.P2
                    switch dim
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

                        otherwise
                            error('MATLAB:ofem:finiteelement:quaddata:InvalidMesh',...
                                  'Only 1D, 2D or 3D meshes supported, but %dD mesh found!',dim);
                    end
                    
                case ofem.finiteelement.Q1
                    error('MATLAB:ofem:finiteelement:weights:NotSupported',...
                          'Q1 not implemented yet!');
                    
                case ofem.finiteelement.Q2
                    error('MATLAB:ofem:finiteelement:weights:NotSupported',...
                          'Q2 not implemented yet!');
                    
                case ofem.finiteelement.NE0P
                    error('MATLAB:ofem:finiteelement:weights:NotSupported',...
                          'NE0P not implemented yet!');

                case ofem.finiteelement.NE1P
                    error('MATLAB:ofem:finiteelement:weights:NotSupported',...
                          'NE1P not implemented yet!');

                case ofem.finiteelement.NE2P
                    error('MATLAB:ofem:finiteelement:weights:NotSupported',...
                          'NE2P not implemented yet!');
            end
        end

        %%
        function pipj = phiiphij(obj,dim)
        %PHIIPHIJ returns the pointwise product of shape functions.
        %
        % pipj=PHIIPHIJ(dim) returns
        % $$
        %   \int_{\Delta} \varphi_i\varphi_j\,dx
        % $$
        % i.e. the integral of the pointwise product of shape functions
        % over the reference entity in the specified dimension dim. The
        % format of pipj is as follows:
        %   pipj(i,j)
        %

            switch obj
                case ofem.finiteelement.P1
                    pipj = (ones(dim+1)+eye(dim+1))/factorial(2+dim);
                    
                case ofem.finiteelement.P2
                    
                    switch dim
                        case 1
                            % node (0,0)        phi1=(2*l1-1)*l1   node 1
                            % node (1,0)        phi2=(2*l2-1)*l2   node 2
                            % node (0,1)        phi3=4*l1*l2       edge 12
                            
                            pipj = [ 4, -1, -1; ...
                                    -1,  4, -1; ...
                                    -1, -1,  4 ]/30;
                            
                        case 2
                            % node (0,0)        phi1=(2*l1-1)*l1   node 1
                            % node (1,0)        phi2=(2*l2-1)*l2   node 2
                            % node (0,1)        phi3=(2*l3-1)*l3   node 3
                            % node (0.5,0)      phi4=4*l1*l2       edge 12
                            % node (0.5,0.5)    phi5=4*l2*l3       edge 23
                            % node (0,0.5)      phi6=4*l1*l3       edge 31
                            
                            pipj = [ 1/60 , -1/360, -1/360,     0, -1/90,     0; ...
                                    -1/360,  1/60 , -1/360,     0,     0, -1/90; ...
                                    -1/360, -1/360,  1/60 , -1/90,     0,     0; ...
                                        0 ,     0 , -1/90 ,  4/45,  2/45,  2/45; ...
                                    -1/90 ,     0 ,     0 ,  2/45,  4/45,  2/45; ...
                                        0 , -1/90 ,     0 ,  2/45,  2/45,  4/45];
                        case 3
                            % node (0,0,0)        phi1 =(2*l1-1)*l1        node 1
                            % node (1,0,0)        phi2 =(2*l2-1)*l2        node 2
                            % node (0,1,0)        phi3 =(2*l3-1)*l3        node 3
                            % node (0,0,1)        phi4 =(2*l4-1)*l4        node 4
                            % node (0.5,0,0)      phi5 =4*l1*l2            edge 12
                            % node (0,0.5,0)      phi6 =4*l1*l3            edge 13
                            % node (0,0,0.5)      phi7 =4*l1*l4            edge 14
                            % node (0.5,0.5,0)    phi8 =4*l2*l3            edge 23
                            % node (0.5,0,0.5)    phi9 =4*l2*l4            edge 24
                            % node (0,0.5,0.5)    phi10=4*l3*l4            edge 34
                            
                            pipj = [ 1/480 , -1/2520, -1/2520, -1/2520, -1/630, -1/630, -1/630, -1/420, -1/420, -1/420; ...
                                    -1/2520,  1/480 , -1/2520, -1/2520, -1/630, -1/420, -1/420, -1/630, -1/630, -1/420; ...
                                    -1/2520, -1/2520,  1/480 , -1/2520, -1/420, -1/630, -1/420, -1/630, -1/420, -1/630; ...
                                    -1/2520, -1/2520, -1/2520,  1/480 , -1/420, -1/420, -1/630, -1/420, -1/630, -1/630; ...
                                    -1/630 , -1/630 , -1/420 , -1/420 ,  4/315,  2/315,  2/315,  2/315,  2/315,  1/315; ...
                                    -1/630 , -1/420 , -1/630 , -1/420 ,  2/315,  4/315,  2/315,  2/315,  1/315,  2/315; ...
                                    -1/630 , -1/420 , -1/420 , -1/630 ,  2/315,  2/315,  4/315,  1/315,  2/315,  2/315; ...
                                    -1/420 , -1/630 , -1/630 , -1/420 ,  2/315,  2/315,  1/315,  4/315,  2/315,  2/315; ...
                                    -1/420 , -1/630 , -1/420 , -1/630 ,  2/315,  1/315,  2/315,  2/315,  4/315,  2/315; ...
                                    -1/420 , -1/420 , -1/630 , -1/630 ,  1/315,  2/315,  2/315,  2/315,  2/315,  4/315 ];
                    end
                    
                    
                case ofem.finiteelement.Q1
                    error('ofem:finiteelement:phiiphij:Unsupported',...
                          'Q1 not implemented yet!');
                    
                case ofem.finiteelement.Q2
                    error('ofem:finiteelement:phiiphij:Unsupported',...
                          'Q2 not implemented yet!');
                    
                case ofem.finiteelement.NE0P
                    error('ofem:finiteelement:phiiphij:Unsupported',...
                          'NE0P not implemented yet!');
                    
                case ofem.finiteelement.NE1P
                    error('ofem:finiteelement:phiiphij:Unsupported',...
                          'NE1P not implemented yet!');
                    
                case ofem.finiteelement.NE2P
                    error('ofem:finiteelement:phiiphij:Unsupported',...
                          'NE2P not implemented yet!');
            end
        end


        %%
        function shape = phi(obj,l)
        %PHI returns the shape functions.
        %
        % shape=PHI(l) returns the shape functions evaluated at the
        % barycentric coordinates given by l. l is expected as returned by
        % quaddata, i.e. l contains per row the barycentric coordinates of
        % the evaluation points.
        %
        % shape has the following format
        %   shape(sidx,lidx)
        % with
        %   sidx: index of the shape function (1 for first, ...)
        %   lidx: index of the evaluation point
        %
        % see also OFEM.FINITEELEMENT.QUADDATA
        %
            switch obj
                case ofem.finiteelement.P1
                    shape = l';
                    
                case ofem.finiteelement.P2
                    dim = size(l,2)-1;
                    
                    switch dim
                        case 1
                            % node (0,0)        phi1=(2*l1-1)*l1   node 1
                            % node (1,0)        phi2=(2*l2-1)*l2   node 2
                            % node (0,1)        phi3=4*l1*l2       edge 12
                            
                            shape = [ (2*l-1).*l, 4*prod(l,2) ]';
                            
                        case 2
                            % node (0,0)        phi1=(2*l1-1)*l1   node 1
                            % node (1,0)        phi2=(2*l2-1)*l2   node 2
                            % node (0,1)        phi3=(2*l3-1)*l3   node 3
                            % node (0.5,0)      phi4=4*l1*l2       edge 12
                            % node (0.5,0.5)    phi5=4*l2*l3       edge 23
                            % node (0,0.5)      phi6=4*l1*l3       edge 31
                            
                            shape = [ (2*l-1).*l,...
                                       4*prod(l(:,[1,2]),2),...
                                       4*prod(l(:,[2,3]),2),...
                                       4*prod(l(:,[3,1]),2) ]';
                        case 3
                            % node (0,0,0)        phi1 =(2*l1-1)*l1        node 1
                            % node (1,0,0)        phi2 =(2*l2-1)*l2        node 2
                            % node (0,1,0)        phi3 =(2*l3-1)*l3        node 3
                            % node (0,0,1)        phi4 =(2*l4-1)*l4        node 4
                            % node (0.5,0,0)      phi5 =4*l1*l2            edge 12
                            % node (0,0.5,0)      phi6 =4*l1*l3            edge 13
                            % node (0,0,0.5)      phi7 =4*l1*l4            edge 14
                            % node (0.5,0.5,0)    phi8 =4*l2*l3            edge 23
                            % node (0.5,0,0.5)    phi9 =4*l2*l4            edge 24
                            % node (0,0.5,0.5)    phi10=4*l3*l4            edge 34
                            
                            shape = [ (2*l-1).*l           ,...
                                       4*prod(l(:,[1,2]),2),...
                                       4*prod(l(:,[1,3]),2),...
                                       4*prod(l(:,[1,4]),2),...
                                       4*prod(l(:,[2,3]),2),...
                                       4*prod(l(:,[2,4]),2),...
                                       4*prod(l(:,[3,4]),2) ]';
                    end
                    
                    
                case ofem.finiteelement.Q1
                    error('ofem:finiteelement:phi:Unsupported',...
                          'Q1 not implemented yet!');
                    
                case ofem.finiteelement.Q2
                    error('ofem:finiteelement:phi:Unsupported',...
                          'Q2 not implemented yet!');
                    
                case ofem.finiteelement.NE0P
                    error('ofem:finiteelement:phi:Unsupported',...
                          'NE0P not implemented yet!');
                    
                case ofem.finiteelement.NE1P
                    error('ofem:finiteelement:phi:Unsupported',...
                          'NE1P not implemented yet!');
                    
                case ofem.finiteelement.NE2P
                    error('ofem:finiteelement:phi:Unsupported',...
                          'NE2P not implemented yet!');
            end
        end


        %%
        function dshape = dphi(obj,l)
        %DPHI returns the gradients of shape functions.
        %
        % dshape=DPHI(l) returns the gradients of the shape functions
        % evaluated at the barycentric coordinate given by l. l is expected
        % as returns by quaddata, i.e. l contains per row the barycentric
        % coordinates of the evaluation points.
        %
        % dshape has the following format
        %   dshape(sidx,didx,lidx)
        % with
        %   sidx: index of the shape function (1 for first, ...)
        %   didx: index of the derivative (1 for x, 2 for y ...)
        %   lidx: index of the evaluation point
        %
        % see also OFEM.FINITEELEMENT.QUADDATA
        %
            
            Nl = size(l,1);
            
            switch obj
                case ofem.finiteelement.P1
                    dim = size(l,2)-1;

                    dshape = repmat([-1*ones(1,dim);eye(dim)],1,1,size(l,1));
                    
                case ofem.finiteelement.P2
                    dim = size(l,2)-1;

                    switch dim
                        case 1
                            % node (0,0): phi1=(2*l1-1)*l1   node 1
                            % node (1,0): phi2=(2*l2-1)*l2   node 2
                            % node (0,1): phi3=4*l1*l2       edge 12

                            l1=l(:,1);
                            l2=l(:,2);

                            dshape = [ -4*l1+1    ; ... %d(phi1)/dx
                                        4*l2-1    ; ... %d(phi2)/dx
                                        4*(l1-l2)];     %d(phi3)/dx

                            dshape = reshape(dshape,3,1,Nl);

                        case 2
                            % node (0,0)    : phi1=(2*l1-1)*l1   node 1
                            % node (1,0)    : phi2=(2*l2-1)*l2   node 2
                            % node (0,1)    : phi3=(2*l3-1)*l3   node 3
                            % node (0.5,0)  : phi4=4*l1*l2       edge 12
                            % node (0.5,0.5): phi5=4*l2*l3       edge 23
                            % node (0,0.5)  : phi6=4*l1*l3       edge 31

                            l1=l(:,1)';
                            l2=l(:,2)';
                            l3=l(:,3)';

                            dshape = [ -4*l1+1      ; ... %d(phi1)/dx
                                        4*l2-1      ; ... %d(phi2)/dx
                                        zeros(1,Nl) ; ... %d(phi3)/dx
                                        4*(l1-l2)   ; ... %d(phi4)/dx
                                        4*l3        ; ... %d(phi5)/dx
                                       -4*l3        ; ... %d(phi6)/dx
                                       -4*l1+1      ; ... %d(phi1)/dy
                                        zeros(1,Nl) ; ... %d(phi2)/dy
                                        4*l3-1      ; ... %d(phi3)/dy
                                       -4*l2        ; ... %d(phi4)/dy
                                        4*l2        ; ... %d(phi5)/dy
                                        4*(l1-l3)  ];     %d(phi6)/dy
                            
                            dshape = reshape(dshape,6,2,Nl);
                            
                        case 3
                            % node (0,0,0)    : phi1=(2*l1-1)*l1   node 1
                            % node (1,0,0)    : phi2=(2*l2-1)*l2   node 2
                            % node (0,1,0)    : phi3=(2*l3-1)*l3   node 3
                            % node (0,0,1)    : phi4=(2*l4-1)*l4   node 4
                            % node (0.5,0,0)  : phi5=4*l1*l2       edge 12
                            % node (0,0.5,0)  : phi6=4*l1*l3       edge 13
                            % node (0,0,0.5)  : phi7=4*l1*l4       edge 14
                            % node (0.5,0.5,0): phi8=4*l2*l3       edge 23
                            % node (0.5,0,0.5): phi9=4*l2*l4       edge 24
                            % node (0,0.5,0.5): phi10=4*l3*l4      edge 34
                            
                            l1=l(:,1)';
                            l2=l(:,2)';
                            l3=l(:,3)';
                            l4=l(:,4)';
                            
                            dshape = [ -4*l1+1     ; ... %d(phi1 )/dx
                                        4*l2-1     ; ... %d(phi2 )/dx
                                        zeros(1,Nl); ... %d(phi3 )/dx
                                        zeros(1,Nl); ... %d(phi4 )/dx
                                        4*(l1-l2)  ; ... %d(phi5 )/dx
                                       -4*l3       ; ... %d(phi6 )/dx
                                       -4*l4       ; ... %d(phi7 )/dx
                                        4*l3       ; ... %d(phi8 )/dx
                                        4*l4       ; ... %d(phi9 )/dx
                                        zeros(1,Nl); ... %d(phi10)/dx
                                       -4*l1+1     ; ... %d(phi1 )/dy
                                        zeros(1,Nl); ... %d(phi2 )/dy
                                        4*l3-1     ; ... %d(phi3 )/dy
                                        zeros(1,Nl); ... %d(phi4 )/dy
                                       -4*l2       ; ... %d(phi5 )/dy
                                        4*(l1-l3)  ; ... %d(phi6 )/dy
                                       -4*l4       ; ... %d(phi7 )/dy
                                        4*l2       ; ... %d(phi8 )/dy
                                        zeros(1,Nl); ... %d(phi9 )/dy
                                        4*l4       ; ... %d(phi10)/dy
                                       -4*l1+1     ; ... %d(phi1 )/dz
                                        zeros(1,Nl); ... %d(phi2 )/dz
                                        zeros(1,Nl); ... %d(phi3 )/dz
                                        4*l4-1     ; ... %d(phi4 )/dz
                                       -4*l2       ; ... %d(phi5 )/dz
                                       -4*l3       ; ... %d(phi6 )/dz
                                        4*(l1-l4)  ; ... %d(phi7 )/dz
                                        zeros(1,Nl); ... %d(phi8 )/dz
                                        4*l2       ; ... %d(phi9 )/dz
                                        4*l3      ];     %d(phi10)/dz
                            
                            dshape = reshape(dshape,10,3,Nl);

                        otherwise
                            error('ofem:finiteelement:dphi:InvalidMesh',...
                                  'Only 1D, 2D or 3D meshes supported, but %dD mesh found!',dim);
                    end
                    
                case ofem.finiteelement.Q1
                    error('ofem:finiteelement:dphi:Unsupported',...
                          'Q1 not implemented yet!');
                    
                case ofem.finiteelement.Q2
                    error('ofem:finiteelement:dphi:Unsupported',...
                          'Q2 not implemented yet!');
                    
                case ofem.finiteelement.NE0P
                    error('ofem:finiteelement:dphi:Unsupported',...
                          'NE0P not implemented yet!');
                    
                case ofem.finiteelement.NE1P
                    error('ofem:finiteelement:dphi:Unsupported',...
                          'NE1P not implemented yet!');
                    
                case ofem.finiteelement.NE2P
                    error('ofem:finiteelement:dphi:Unsupported',...
                          'NE2P not implemented yet!');
            end
        end
    end
end

