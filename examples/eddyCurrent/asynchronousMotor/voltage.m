function Y=voltage(X,coil,Utot)

Y=0*X(1,:,:);

for i=1:size(coil.radius,1)

    XX = double(permute(X,[3,1,2]));
    XX(:,1) = XX(:,1)-coil.pos(i,1);
    XX(:,2) = XX(:,2)-coil.pos(i,2);
    idx = double(dot(XX,XX,2))<=coil.radius(i)^2;
%     XX(:,1)=X(:,1)-coil.pos(i  ,1);
%     XX(:,2)=X(:,2)-coil.pos(i  ,2);
%     idx  = dot(XX,XX,2)<=coil.radius(i  )^2;

    if ~isempty(idx)
        Y(1,:,idx) = Utot;
    end

end

end