

classdef nedelecP < handle
	properties
		order
		N
		curlN
	end
	
	methods
		
		%%
		function obj = nedelecP(order)
			obj.order = order;
		end
		
		%%
		function [w,l] = quaddata(~,dim,order)
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
					l = [l1,l2,l3,l4];
					
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
					
					N0 = [N1,N2,N3,N4,N5,N6];
					
					switch obj.order
						case 1
							%Order 1 just uses the Ns
							baseH = N0;
						case 2
							% Order 2 uses Edge and Face elements. Code is
							% ugly but fast
							idxE = repelem(1:6,1,2);
							idxB = repelem([1,2,1,3,1,4,2,3,2,4,3,4],3,1);
							baseE = N0(:,idxE).*l(idxB);
							idxF = [1,4,1,5,4,6,2,6];
							idxB = repelem([3,1,4,1,4,2,4,1],3,1);
							baseF = N0(:,idxF).*l(idxB);
							baseH = [baseE,baseF];
						otherwise
							% For all higher orders we need edge, face and
							% interior functions. For now this is a big
							% TODO
							%% Start with the edges
							for i = 1:6
								for j = 1:obj.order
								end
							end
					end
					
					%% Next we need to check the required order and compute the corresponding basis
					%% Fix this from here!
					% 				H1 = l1*N1;
					% 				H2 = l2*N1;
					% 				H3 = l1*N2;
					% 				H4 = l3*N2;
					% 				H5 = l1*N3;
					% 				H6 = l4*N3;
					% 				H7 = l2*N4;
					% 				H8 = l3*N4;
					% 				H9 = l2*N5;
					% 				H10 = l4*N5;
					% 				H11 = l3*N6;
					% 				H12 = l4*N6;
					% 				H13 = l3*N1;
					% 				H14 = l1*N4;
					% 				H15 = l4*N1;
					% 				H16 = l1*N5;
					% 				H17 = l4*N4;
					% 				H18 = l2*N6;
					% 				H19 = l4*N2;
					% 				H20 = l1*N6;
					%
					% 				baseH = [H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,H15,H16,H17,H18,H19,H20];
					% 				basecH = [cH1,cH2,cH3,cH4,cH5,cH6,cH7,cH8,cH9,cH10,cH11,cH12,cH13,cH14,cH15,cH16,cH17,cH18,cH19,cH20];
			end
			% Compute the curl of the generated basis functions
			basecH = baseH;
			for i = 1:size(baseH,2)
				basecH(:,i) = curl(baseH(:,i),dr);
			end
			[w,l] = obj.quaddata(3,3);
			obj.N = repmat(baseH,1,1,size(l,1));
			obj.curlN = repmat(basecH,1,1,size(l,1));
			for i=1:size(obj.N,3)
				obj.N(:,:,i) = (subs(obj.N(:,:,i),dr,l(i,:)));
				obj.curlN(:,:,i) = (subs(obj.curlN(:,:,i),dr,l(i,:)));
			end
			obj.N = double(obj.N);
			obj.curlN = double(obj.curlN);
			N = obj.N;
			curlN = obj.curlN;
			
		end
	end
end














