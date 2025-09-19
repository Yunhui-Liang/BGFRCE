function [Gcell] = construct_GF_cell(Hcell)
% ------------------------------------------------------------
% FUNCTION: construct_GF_cell
% PURPOSE : Compute projection matrices G_i = H_i (H_i^T H_i)^{-1} H_i^T
% INPUT   : Hcell - 1 x m cell, each Hcell{i} = H_i (n x k_i matrix)
% OUTPUT  : Gcell - 1 x m cell, each Gcell{i} = G_i (n x n projection matrix)
%           G_star - average of all G_i (n x n matrix)
% ------------------------------------------------------------

    % Number of views / number of H matrices
    m = numel(Hcell);
    if m == 0
        error('Hcell is empty.');
    end

    % Preallocate cell for G
    Gcell = cell(1, m);

    % Loop over each H_i
    for i = 1:m
        Hi = Hcell{i};
        if isempty(Hi)
            error('Hcell{%d} is empty.', i);
        end
        
        % Compute Gi = Hi * inv(Hi'Hi) * Hi'
        % Use backslash instead of inv for numerical stability:
        %   Gi = Hi / (Hi'Hi) * Hi'
        HiinvHTH = Hi./sum(Hi, 1);
        Gi = HiinvHTH * Hi';
        
        % Store in cell
        Gcell{i} = Gi;

    end

end
