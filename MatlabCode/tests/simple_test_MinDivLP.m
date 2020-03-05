%% A very simple test showing how to test the MinDivLP code

% add the code repository to the Matlab path
addpath(strcat(fileparts(pwd), '/src'))

% create some simple data
A_k_small = [.5, 0, 0; 
            0, 1/3, 1/5;
            .5, 2/3, 4/5];

A_k_large = [0 ,0 ,1/6;
            1/4, 1/3, 0;
            1/4, 0 , 2/6;
            1/4, 1/3, 2/6;
            1/4, 1/3, 1/6];

x_true = [1;0;0];  % to test against

y_small = A_k_small*x_true;
y_large = A_k_large*x_true;

lambda = 10000;
q = 0.1;

% run the algorithm
x_star = MinDivLP(A_k_small, A_k_large, y_small, y_large, lambda,q);


% Then do the tests

% a realistic test
thresh = .00001;
L1_error = norm(x_star - x_true,1);
assert(L1_error < thresh, sprintf('L1 reconstruction error is above %f', thresh))

% a test for EXACT recovery (which will probably never pass on realistic data)
assert(all(x_star == x_true), 'Reconstructed and true vectors are not exactly the same')
