function [ alpha, u ] = modified_gradient_projection(X, Y, C, k, kernel, tol, ...
                                                  maxiter)
    % ---------------------------------------------------------------
    
    % MODIFIED GRADIENT PROJECTION
    %
    % Description:
    %
    %
    %
    % Contributors:
    %
    % - José Manuel Proudinat Silva
    % - Joaquín Sánchez García
    %
    %
    % Input:
    %
    % - X : [ data matrix : (n, m) ]
    % - Y : [ labels matrix : (n, 1) ]
    % - kernel : [ string : kind of kernel ]
    % - tol : [ float : convergency tolerance ]
    % - maxiter : [ int : maximum number of iterations ]
    %
    %
    % Output:
    %
    % - alpha : [ matrix (n, 1) : dual problem result ]
    %
    %
    % Modelos Matemáticos y Numéricos (ITAM)
    % 2015
    %
    % ---------------------------------------------------------------

    % Start with an empty working set
    n = length(Y);
    % Initial value alpha = 0
    alpha = sparse(n, 1);
    % Empty kernel matrix
    K = sparse(n, n);
    K = X' * X;
    u = alpha;
    g = ones(n, 1);
    rho = mean(Y);   % Heurística Proudinat-Sánchez 1
    iter = 1;
    
    % While there are examples violating the optimzality condition
    tic;

    % TODO : Elimiar esto cuando tenga la función que agrega vectores al
    % working set
    [optim, B] = ChecaOptim(alpha, Y, g, rho, C, k)
    
    while not(optim) && iter < maxiter
  
        %
        % loop
        %
        loop = 1;
        while loop
            % Compute the gradient
            aux = Y .* alpha;
            aux = K * aux;
            aux = Y .* aux;
            g(B) = 1 - aux(B)

            %
            % repeat
            %

            %
            % Projection enviction loop
            %
            change = 1;
            while change
                % Initialize the direction
                % TODO: Esto debe optimizarse para no calcular nuevamente 
                % las derivadas parciales que ya estaban
                u = sparse(n, 1);
                
                % Compute the rho from the optimality condition as a mean
                aux = Y(B) .* g(B);
                rho = mean(aux)
        
                % Compute the next direction
                u(B) = g(B) - rho * Y(B)       
                
                % Search for positive and negative directions and second condition
                B2 = find(u > 0 & (C - alpha) < tol)
                B3 = find(u < 0 & alpha < tol)
                index = union(B2, B3)
                if isempty(index)
                    aux = B;
                else
                    aux = setdiff(B, index);
                end
                if isempty(aux)
                    change = 0;
                    u = sparse(n, 1);
                end
                
                if length(aux) == length(B)
                    change = 0;
                end
                B = sort(aux);
            end

            if norm(u, inf) < tol
                fprintf('Ya no hay ascenso lineal \n')
                loop = 0;
            else
                fprintf('Sí hubo \n')
            
                %
                % Compute lambda*
                %
                
                pos = find(u > 0);
                neg = find(u < 0);    

                if not(isempty(pos))
                    lambda_max = min((C - alpha(pos)) ./ u(pos));
                else
                    lambda_max = inf;
                end               
                lambda_max
                if not(isempty(neg))
                    lambda_max = min(min(-alpha(neg) ./ u(neg)), lambda_max);               
                end
                if isempty(lambda_max)
                    lambda_max = 0;
                end
                
                % LLamamos al programa que calcula con Newton
                lambda = LambdaDeNewton(u, g, Y, K)
                % Obtenemos la más pequeña de éstas
                lambda = min(lambda_max, lambda);
                % Obtenemos el tamaño final del paso
                lambda = max(0, lambda)
                
                if not(lambda)               
                    break
                end          

                % Update
                alpha = alpha + lambda * u
                % K = MatrizDeKernel(K, X, 1, alpha);

            end

            iter = iter + 1;

        end

        % Check optimality conditions
        % TODO : Probablemente cuando esté lista la función para elegir el
        % working set no será necesario calcular B en este paso
        [optim, B] = ChecaOptim(alpha, Y, g, rho, C, k)
        if length(B) == 1
            alpha(B) = 1 / Y(B) * (Y' * alpha - Y(B) * alpha(B));
            optim = 1;
        end
        alpha
    end
    
end

        
        
        