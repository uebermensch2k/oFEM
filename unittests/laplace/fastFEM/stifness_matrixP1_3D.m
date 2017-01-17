function [A,areas]=stifness_matrixP1_3D(elements,coordinates)
NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!  
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB           
        coord(d,i,:)=coordinates(elements(:,i),d);       
    end   
end
clear coordinates 
IP=[1/4 1/4 1/4]';    
[dphi,jac] = phider(coord,IP,'P1'); %integration rule, it must be known!  
clear coord 

areas=abs(squeeze(jac))/factorial(DIM);
clear jac

dphi = squeeze(dphi); 
Z=astam(areas',amtam(dphi,dphi));
clear dphi

Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);
clear elements
%copy this part for a creation of a new element

X=permute(Y,[2 1 3]); 
A=sparse(X(:),Y(:),Z(:)); 
