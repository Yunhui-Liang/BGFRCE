function [M2, M] = compute_graph_power_eig_v3(S, order)
% -------------------------------------------------------
% FUNCTION: compute_graph_power_eig
% PURPOSE : Efficiently compute sqrt(M), where
%           M = S + S^2 + ... + S^order
% METHOD  : Eigen-decomposition + element-wise power on eigenvalues
%
% INPUT  :
%   S     - (n x n) symmetric adjacency matrix
%   order - highest order to compute (e.g., 3)
% OUTPUT :
%   M_sqrt - sqrt(M), used for high-order graph filtering
%
% EXAMPLE:
%   M_sqrt = compute_graph_power_eig(S, 3);
%
% -------------------------------------------------------

    if nargin < 2
        order = 5;  % Default to 3rd order
    end

    % === Step 1. Eigen-decomposition ===
    [V, D] = eig(full(S));  % V: eigenvectors, D: diag(eigenvalues)
    eigvals = diag(D);
    % 避免后续计算出现复数
    eigvals = max(eigvals, 0);
    % === Step 2. Efficiently compute sum of powers ===
    eig_sum = zeros(size(eigvals));

    % Use vectorized element-wise power
    for p = 1:order
        eig_sum = eig_sum + eigvals.^p;
    end

    % === Step 3. Compute M^2 in eigen-space ===
    eig_sum2 = eig_sum.^2;

    % === Step 4. Reconstruct M ===
    M = V * diag(eig_sum) * V';
    M = (M + M')/2;
    M2 = V * diag(eig_sum2) * V';
    M2 = (M2 + M2')/2;
end
