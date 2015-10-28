function [ optim,B ] = ChecaOptim( alfa,Y,g,rho,C,k,tol )


% Este programa checa si alfa es optimo (considerando rho)

% Creator : Joaquín Sánchez García
% Modified by: José Manuel Proudinat Silva

index1 = find( g-rho*Y > 0);
index2 = find( g-rho*Y < 0);

a = C*ones(length(index1),1)-alfa(index1);
b = alfa(index2);


B = [];

if norm( a,Inf) < tol && norm(b,Inf)< tol 
    optim = 1;
    return
else
    
    
    ind1 = a > tol;
    ind2 = b > tol;
    index1 = index1(ind1);
    index2 = index2(ind2); 
    B = sort([index1;index2]);
    optim = 0;
end



k = min(length(B), k);
% We search for the coordinates for which the optimality condition is most
% violated
%aux = [a(index1) , b(index2)];
%aux = sort(aux);
%aux = aux(1:k);
%[~,ind1] = ismember(a,aux);
%[~,ind2] = ismember(b,aux);
%ind1 (ind1 == 0) = [];
%ind2 (ind2 == 0) = [];
%index1 = index1(ind1);
%index2 = index2(ind2);
%B = sort([index1;index2])

 B = B(1:k);

end

