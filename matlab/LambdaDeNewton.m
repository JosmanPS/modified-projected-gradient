function [ lambda ] = LambdaDeNewton(u,g,Y,K)
% Este programa usa la formula de Newton para encontrar 
% el tamaño de paso optimo en la busqueda lineal en la caja

% Al usar la formula de Newton aparece la Hessiana H
% Hij = yi yj Kij donde K es la matriz de Kernel
% pero solo para el producto u'*H*u por lo que solo hay que 
% calcular los elementos del kernel que correspondan a valores
% de u no cero.

ug = u' * g;
Y = sparse(diag(Y));
H = Y * K * Y;
uHu = u' * H * u;

lambda = ug / uHu;

end

