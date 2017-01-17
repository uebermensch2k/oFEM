function M=mass_matrixP1_2D(elements,areas)

Xscalar=kron(ones(1,3),elements); Yscalar=kron(elements,ones(1,3)); 
Zmassmatrix=kron(areas,reshape((ones(3)+eye(3))/12,1,9));
M=sparse(Xscalar,Yscalar,Zmassmatrix);