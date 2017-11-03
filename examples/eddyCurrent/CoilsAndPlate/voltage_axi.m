function Y=voltage_axi(X,coil,Utot,R)

Y=0*X(1,:,:);

Rtot=0;

for i=1:size(coil.min,1)
    idx = coil.min(i,1)<=X(1,:,:) & coil.max(i,1)>=X(1,:,:) & ...
          coil.min(i,2)<=X(2,:,:) & coil.max(i,2)>=X(2,:,:);
      idx = squeeze(idx);

    if ~isempty(idx)
%         Y(idx) = Utot/size(coil.min,1);
        h      = coil.max(i,2)-coil.min(i,2);
        Rk     = 1/(h*log(coil.max(i,1)/coil.min(i,1)));
        rho    = X(1,:,idx);
        
        if ~isempty(rho)

          Y(:,:,idx) = Utot*Rk./(2*pi*rho);
        end

        Rtot   = Rtot+Rk;
    end
end

if Rtot ~= 0
    Y=Y./Rtot;
end
end