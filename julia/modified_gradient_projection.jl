function modified_projected_gradient(X::Array{Float64}, 
									 Y::Array{Int64}; 
									 C::Float64=1.0, 
									 kernel::ASCIIString="linear", 
									 tol::Float64=1e-4, 
									 maxiter::Int64=100)

	#=

	          MODIFIED PROJECTED GRADIENT

	Description
	------------

	Contributors
	-------------
	- José Manuel Proudinat Silva
	    jmps2812@gmail.com
	- Joaquín Sánchez García
	    joaqsangar17@gmail.com

	Input
	------
	-
	-

	Output
	-------
	-
	-

	=#

	# Initial values
	# ---------------
	# TODO : Add starting time variable 'tic'
    n = length(Y)
    alpha = spzeros(n, 1)      # Initial value alpha = 0
    K = X' * X                 # TODO : Compute only necessary parts of K
    u = spzeros(n, 1)
    g = ones(n, 1)
    rho = mean(Y)
    iter = 1

    # While there are examples violating th optimality condition...
    optim, B = CheckOptim(alpha, Y, g, rho, C, k)

    while !optim && iter < maxiter

    	loop = true
    	while loop

    		# Compute the gradient
    		# ---------------------
    		aux = Y .* alpha
    		aux = K * alpha
    		aux = Y .* aux
    		g[B] = 1 - aux[B]

    		# Projection enviction loop
    		# --------------------------
    		change = true
    		while change

    			# Initialize the direction
    			u = zeros(n, 1)

    			# Compute the rho from the optimality condition
    			# as a mean
    			aux = Y[B] .* g[B]
    			rho = mean(aux)

    			# Compute the next direction
    			u[B] = g[B] - rho * Y[B]

    			# Search for positive and negative directions 
    			# and second condition
    			

end
