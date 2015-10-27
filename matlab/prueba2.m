X = [5 1 6 1.5; 7 2 8 3; 2 8 1.5 7.5]/10;
% X = X';
Y = [1; -1; 1; -1];

modified_gradient_projection(X, Y, 1, 3, 1, 1e-4, 200)