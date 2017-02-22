function Y = current(X)

rho = sqrt(dot(X,X,2)); 
I   = 1e-4;
Q   = pi*0.001^2;

Y   = zeros(size(X,1),1);
Y(rho<0.001) = I/Q;

end
