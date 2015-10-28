function [ optim,B ] = ChecaOptim( alfa,Y,g,rho,C,k,tol )


% Este programa checa si alfa es optimo (considerando rho)

% Creator : Joaquín Sánchez García
% Modified by: José Manuel Proudinat Silva

index1 = find( g-rho*Y > 0 );
index2 = find( g-rho*Y < 0 );

a = C - alfa(index1);
b = alfa(index2);


B = [];

if norm(a, Inf) < tol && norm(b, Inf)< tol 
    optim = 1;
    return
else    
    ind1 = a > tol;
    ind2 = b > tol;
    a = a(ind1);
    b = b(ind2);
    index1 = index1(ind1);
    index2 = index2(ind2); 
    B = [index1; index2];
    optim = 0;
end

% Get the elements that most violated the optimality condition
violated = [a; b];
[~, order_B] = sort(violated, 'descend');
B = B(order_B);

% Get just k elements
k = min(length(B), k);
B = B(1:k);

end

