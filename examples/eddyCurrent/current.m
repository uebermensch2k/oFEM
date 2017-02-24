function Y = current(X)

rho = squeeze(sqrt(dot(X,X,1)));
I   = 1e-4;
Q   = pi*0.001^2;

Y = 0*X(1,:,:);
Y(:,:,rho<0.001) = I/Q;

end
