function t=tr(A, dim)
%TR returns the trace of a tensor
t=sum(A(:,1:dim+1:end,:),2);
end
