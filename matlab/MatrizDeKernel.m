function [ K ] = MatrizDeKernel( K, X, kernel, alfa,tol)
% Este programa calcula la matriz de Kernel segun el tipo
% aprovecha los valores para los que alfa = 0 y esos no los calcula
% porque siempre que se necesita Kij va multiplicado por alfaj

% Creator : Joaquín Sánchez García
% Modified by: José Manuel Proudinat Silva

index = find (alfa > tol );

if kernel == 1
   
    for i = index
        for j = index
            % Only updates zero values
            if K(i,j) == 0
                K(i,j) = X(:,i)'*X(:,j);
            end            
        end
    end
else
    error('todavia no programo ese kernel');
end


end

