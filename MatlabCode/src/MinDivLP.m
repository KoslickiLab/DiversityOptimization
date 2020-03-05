function x_star = MinDivLP(A_k_small, A_k_large, y_small, y_large, lambda, q)
%MinDivLP
% A basic, regularized version of the MinDivLP algorithm.
% Call via:
% x_star = MinDivLP(A_k_small, A_k_large, y_small, y_large, lambda, q)
%
% Parameters are:
%   A_k_small is the [m_small, N]-sized sensing matrix
%   A_k_large is the [m_large, N]-sized sensing matrix
%   y_small is the data vector of size [m_small, 1]
%   y_large is the data vector of size [m_large, 1]
%   lambda is the regularization paramater (larger values indicated better
%     fit to constraints, at the cost potentially higher execution time and 
%     may lead to over-fitting if set too large. Typical value is 10000 or
%     1000
%   q is the parameter used in the MinDivLP algorithm. Must have 0<q<1,
%     typically, q is set to something like q = 0.1
%
% Returns:
% x_star: an [N, 1] vector

    % Generate the B matrix
    B = (A_k_large > 0);  % create the B matrix
    f = 1./(B'*y_large).^(1-q);  % weighting factor in MinDivLP, no noise in this y vector

    % Do the Quikr trick in MinDivLP and bring the A^{(h)} y^{(h)} = x inside the
    % objective function via lambda^2*||A^{(h)} y^{(h)} - x||_2^2
    % then we can use the active set algorithm lsqnonneg
    x_star = lsqnonneg([f'; lambda*A_k_small], [0;lambda*y_small]);  % the regularized version of MinDivLP for noisy data
    x_star = x_star./sum(x_star);  % normalize
end