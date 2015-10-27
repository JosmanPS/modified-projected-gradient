function [ optim,B ] = ChecaOptim( alfa,Y,g,rho,C,k )


% Este programa checa si alfa es optimo (considerando rho)

% Creator : Joaquín Sánchez García
% Modified by: José Manuel Proudinat Silva

index1 = find( g-rho*Y > 0);
index2 = find( g-rho*Y < 0);

a = C*ones(length(index1),1)-alfa(index1);
b = alfa(index2);

viola = [index1' index2'; a' b']';

B = [];

if norm( a,Inf) < 1.e-8 && norm(b,Inf)<1.e-8
    optim = 1;
    return
else
    ind1 = find(a > 1.e-8);
    ind2 = find(b > 1.e-8);
    index1 = index1(ind1);
    index2 = index2(ind2);
    B = sort([index1;index2]);
    optim = 0;
end

k = min(length(B), k);
B = B(1:k);

end

