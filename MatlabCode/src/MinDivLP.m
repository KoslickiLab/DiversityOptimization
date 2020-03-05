function x_star = MinDivLP(A_k_small, A_k_large, y_small, y_large, lambda, q)
    % Generate the B matrix
    B = (A_k_large > 0);  % create the B matrix
    f = 1./(B'*y_large).^(1-q);  % weighting factor in MinDivLP, no noise in this y vector

    % Do the Quikr trick in MinDivLP and bring the A^{(h)} y^{(h)} = x inside the
    % objective function via lambda^2*||A^{(h)} y^{(h)} - x||_2^2
    % then we can use the active set algorithm lsqnonneg
    x_star = lsqnonneg([f'; lambda*A_k_small], [0;lambda*y_small]);  % the regularized version of MinDivLP for noisy data
    x_star = x_star./sum(x_star);  % normalize
end

