

classdef nedelecP
    properties
		order
		N
		curlN
	end

    methods
		
		%%
		function obj = nedelecP(order,dim)
			obj.order = order;
		end
		
        %%
        function [w,l] = quaddata(obj,dim,order)
        %QUADDATA returns quadrature points and weights.
        %
        % [w,l]=QUADDATA(dim) returns per row the barycentric coordinates l
        % of the quadrature points and weights w sufficient for the
        % specified element and space dimension given through dim. The
        % quadrature rule is a Gaussian quadrature. The quadarture points l
        % yield barycentric coordinates as they are used by
        % phi and dphi. The weight already contains the jacobian scaling
        % factor for the length/area/volume
        %
        % see also OFEM.NEDELECP.BASIS
        %
			switch dim
				case 2
					w=1/factorial(dim);
					l=repmat(1/(dim+1),1,dim+1);
					coord = [0 1 0; 0 0 1];
					l = coord*l';
				case 3
					switch order
						case 1
							w = 1/6;
							l = [1/4,1/4,1/4];
						case 2
							w=[1;1;1;1]/24;
							a = 0.1381966011250105151795413165634361;
							b = [a,a,a,1-3*a;a,a,1-3*a,a;a,1-3*a,a,a;1-3*a,a,a,a];
							coord = [0 1 0 0; 0 0 1 0; 0 0 0 1];
							l(:,1) = coord*b(:,1);
							l(:,2) = coord*b(:,2);
							l(:,3) = coord*b(:,3);
							l(:,4) = coord*b(:,4);
						case 3
							w = [-4/5,9/20,9/20,9/20,9/20]/6;
							l = [1/4,1/4,1/4;
								1/6,1/6,1/6;
								1/2,1/6,1/6;
								1/6,1/2,1/6;
								1/6,1/6,1/2];
					end
			end
		end


        %%
        function [N,curlN] = basis(obj,dim)
        %BASIS returns the nedelec elements and their curl evaluated at the
        %corresponding quadrature points. Here N are the basis functions
        %evaluated at the quadrature points for order p and curlN is the
        %curl of the basis functions evaluated at quadrature points exact
        %to p-1.
		
		switch dim
			case 2
				error('ofem:finiteelement:basis:Unsupported',...
                          '2D elements are not implemented yet!');
			case 3
				%% First compute the zeroth order Nedelec elements to construct everything
				syms u v w;
				dr = [u,v,w];

				%% H^1 Basis
				l1 = 1-u-v-w;
				l2 = u;
				l3 = v;
				l4 = w;

				%% Gradient of H^1 Basis
				dl1 = gradient(l1,dr);
				dl2 = gradient(l2,dr);
				dl3 = gradient(l3,dr);
				dl4 = gradient(l4,dr);

				%% Zeroth order Nedelec
				N1 = l1*dl2-l2*dl1;
				N2 = l1*dl3-l3*dl1;
				N3 = l1*dl4-l4*dl1;
				N4 = l2*dl3-l3*dl2;
				N5 = l2*dl4-l4*dl2;
				N6 = l3*dl4-l4*dl3;

				%% Next we need to check the required order and compute the corresponding basis
				%% Fix this from here!
				H1 = l1*N1;
				H2 = l2*N1;
				H3 = l1*N2;
				H4 = l3*N2;
				H5 = l1*N3;
				H6 = l4*N3;
				H7 = l2*N4;
				H8 = l3*N4;
				H9 = l2*N5;
				H10 = l4*N5;
				H11 = l3*N6;
				H12 = l4*N6;
				H13 = l3*N1;
				H14 = l1*N4;
				H15 = l4*N1;
				H16 = l1*N5;
				H17 = l4*N4;
				H18 = l2*N6;
				H19 = l4*N2;
				H20 = l1*N6;

				%% Curl of first order Nedelec
				cH1 = curl(H1,dr);
				cH2 = curl(H2,dr);
				cH3 = curl(H3,dr);
				cH4 = curl(H4,dr);
				cH5 = curl(H5,dr);
				cH6 = curl(H6,dr);
				cH7 = curl(H7,dr);
				cH8 = curl(H8,dr);
				cH9 = curl(H9,dr);
				cH10 = curl(H10,dr);
				cH11 = curl(H11,dr);
				cH12 = curl(H12,dr);
				cH13 = curl(H13,dr);
				cH14 = curl(H14,dr);
				cH15 = curl(H15,dr);
				cH16 = curl(H16,dr);
				cH17 = curl(H17,dr);
				cH18 = curl(H18,dr);
				cH19 = curl(H19,dr);
				cH20 = curl(H20,dr);

				baseH = [H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,H15,H16,H17,H18,H19,H20];
				basecH = [cH1,cH2,cH3,cH4,cH5,cH6,cH7,cH8,cH9,cH10,cH11,cH12,cH13,cH14,cH15,cH16,cH17,cH18,cH19,cH20];
		end
		
		end
    end
end














