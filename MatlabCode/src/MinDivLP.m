%% import the data
% set variables
small_k = 4;  % smaller k-mer size
large_k = 6;  % larger k-mer size

A_k = load(sprintf('../data/97_otus_subset.fasta_A_%d.mat', large_k));
A_k_large = A_k.A_k;

A_k = load(sprintf('../data/97_otus_subset.fasta_A_%d.mat', small_k));
A_k_small = A_k.A_k;

clear('A_k')  % get rid of unneeded variable

%% sub-select the data so things run quickly

cols_vs_rows = 3;  % fix 3-times more columns than rows
num_species = cols_vs_rows*4^small_k;  % reduce the number of columns of 
% the sensing matrix so pictures will be generated in a reasonable amount
% of time.
A_k_small = A_k_small(:, 1:num_species);  % Note: these are already column-normalized to be 1
A_k_large = A_k_large(:, 1:num_species);

%% Generate the B matrix
B = (A_k_large > 0);

%% Set some variables
q = .1;  % fixed, small q value s.t. 0<q<1
support_size = 15;  % number of non-zero entries in the simulated ground truth

%% create the simulated ground truth
% create
supp = datasample(1:num_species, support_size, 'Replace', false);  % location of the support
true_x = zeros(num_species,1);  % the true x vector we are trying to reconstruct
true_x(supp) = rand(support_size,1);  % populate with random data
true_x = true_x./sum(true_x);  % normalize to be a probability vector

% noisless y-vectors
y_small_true = A_k_small*true_x;
y_large_true = A_k_large*true_x;

% noisy y-vectors
slop_factor = 2;  % this is for the non-regularized MinDivLP on noisy data
noise_eps = .00001;  % size of noise to add
y_small_noise = A_k_small*true_x + noise_eps*abs(rand(size(A_k_small,1),1));  % add only noise to the small y-vector
y_small_noise = y_small_noise./sum(y_small_noise);

% optimization parameters
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off','Diagnostics','off', 'ConstraintTolerance', .0000001, 'OptimalityTolerance', .0000001);

%% Noisless computations

f = 1./(B'*y_large_true).^(1-q);  % weighting factor in MinDivLP

% MinDivLP

% Using matlab built in linprog: too slow
%[x_star, ~, ~, ~, ~] = linprog(f, [], [], A_k_small, y_sm, zeros(1,num_species), ones(1,num_species), options);  % if you do not have a Gurobi license, you can use the built in linprob
%x_star = x_star/sum(x_star);  % normalize

% Using the linprog_gurobi (code optained from Gurobi). Requires a Gurobi
% license from https://www.gurobi.com/downloads/free-academic-license/ and
% downloading their software from https://www.gurobi.com/downloads/gurobi-optimizer-eula/
[x_star, ~, ~, ~, ~] = linprog_gurobi(f, [], [], A_k_small, y_small_true, zeros(1,num_species), ones(1, num_species), options);  % again, can substitute linprog
x_star = x_star/sum(x_star);  % normalize

%% Noisy computations

f = 1./(B'*y_large_true).^(1-q);  % weighting factor in MinDivLP, no noise in this y vector

% Relaxing the A^{(h)} y^{(h)} = x
% Turns out this is not so good
%[x_star, ~, ~, ~, ~] = linprog_gurobi(f, [A_k_small; -A_k_small], [y_small_noise+slop_factor*noise_eps; -y_small_noise+slop_factor*noise_eps], [], [], zeros(1,num_species), ones(1, num_species), options);  % again, can substitute linprog
%x_star = x_star./sum(x_star);  % normalize

% Instead, do the Quikr trick and bring the A^{(h)} y^{(h)} = x inside the
% objective function via lambda^2*||A^{(h)} y^{(h)} - x||_2^2
% then we can use the active set algorithm lsqnonneg
lambda = 10000;  % technically, this is lambda^2, but whatever
x_star = lsqnonneg([f'; lambda*A_k_small], [0;lambda*y_small_noise]);  % the regularized version of MinDivLP for noisy data
x_star = x_star./sum(x_star);  % normalize
%% Measure the reconstruction accuracy
error_l1 = norm(x_star - true_x, 1);
error_l2 = norm(x_star - true_x, 2);
fprintf('L1 error is: %f\n', error_l1);
fprintf('L2 error is: %f\n', error_l2);

%% Plot things just to sanity check
f = figure();
hold on
plot(1:num_species, x_star, 'bo-');
plot(1:num_species, true_x, 'ko-');
legend('Reconstructed', 'True')

