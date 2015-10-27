function [ alfa ] = MPG(X, Y, C, kernel)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------ ITAM - Modelos Matematicos----------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------ Modified Projected Gradient (Dual Formulation)
% --------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Creator : Joaquín Sánchez García
% Modified by: José Manuel Proudinat Silva



% Iniciamos con valores iniciales. 
% Necesitamos alfa en la caja de tamaño C

alfa = zeros(length(Y),1); 

% alfa=0 no puede ser optima si existe yi <0 .
% el caso en el que no exista no tiene sentido. 
% (No hay nada que clasificar)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Empezamos con working set todos los valores

B = 1:length(alfa);

% El primer valor de los gradientes es puros unos 

g = ones(length(alfa),1); 

% Como todos los valores de alfa son cero la primera matriz de kernel no
% afecta
K = zeros(length(alfa),length(alfa));



% Vamos a usar la variable optim para checar optimalidad en otro programa
% cuando alfa sea optima optim será 1 
% cuando alfa no sea optima optim sera -1

optim = 0; 
tic;
% Entr al loop siempre que alfa no sea optima
while optim == 0 
    
    
    % Empezamos poniendo como ceros en la direccion de avance

    u = zeros(length(alfa),1);
    
    
    
    % Construimos los pasos para escribir el gradiente
    aux = Y.*alfa;
    aux = K*aux;
    aux = Y.*aux;
    
    % Calculamos los gradientes
    g(B) = ones(length(B),1) - aux(B);
    
    aux = Y(B).*g(B);
    
    % Calculamos rho de la condicion de optimalidad como la media
    rho = mean(aux);
    % Calculamos la siguiente dirección 
    u(B) = g(B) - rho*Y(B);
    
    
    if u == 0 
        fprintf ('Ya no hay hay ascenso lineal, el valor de alfa es %3.2e \n',alfa);
    end
    
    % Buscamos cuales de los indices en B tienen u>0 o u<0
    B2 = find( u > 0);
    B3 = find( u < 0);
    % De esos indices checamos cuales cumplen la segunda condicion
    B2alfa = alfa(B2) == C;
    B3alfa = alfa(B3) == 0;
    % Recuperamos los indices originales
    B2 = B2(B2alfa);
    B3 = B3(B3alfa);
    % Acomodamos estos dos vectores 
    vec = [B2;B3];
    % Reacomodamos vec 
    vec = sort(vec);
    % Identificamos en B las componentes que cumplieron alguna de las dos
    % condiciones poniendole 0 a su entrada
    B(vec) = 0; 
    % Quitamos de B las entradas correspondientes
    B (B == 0) = [];
    
        
    % Calculamos el maximo valor del paso 
    index = find(u > 0);
    index2= find(u < 0);
    
    if isempty(index2) == 1
        lambdamax = min( (C*ones(length(index),1) - alfa(index))./(u(index)));
    else
      lambdamax = min(min(-alfa(index2)./u(index2)),min( (C*ones(length(index),1) - alfa(index))./(u(index))));
    end
    
    if isempty(lambdamax) ==1
          lambdamax =0;
          
    end
      
    % Calculamos el mayor valor de lambda para que no se salga de la caja
   
    % Llamamos al programa que calcula con Newton.
    [lambdaD] = LambdaDeNewton(u,g,Y,K);
   
    
    % Obtenemos la mas pequeña de estas
    minlambda = min(lambdamax,lambdaD);
    
    
    % Obtenemos el tamaño final del paso
    lambda = max(0,minlambda);
    
   if lambda == 0
        toc
       return;
      
   end
       
    
    alfa = alfa+lambda*u;
    
    [K] = MatrizDeKernel(K, X,1,alfa);
    
    % Programa que checa el criterio de optimalidad alfa rho
      
    [optim , B ] = ChecaOptim(alfa,Y,g,rho,C, 3);
    
    
    
end


end



