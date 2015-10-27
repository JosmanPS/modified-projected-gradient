X = csvread('data_test.csv', 1, 1);

Y = X(:, 3);
X = X(:, 1:2);
X = X';


modified_gradient_projection(X, Y, 1000, 50, 1, 1e-4, 2000)