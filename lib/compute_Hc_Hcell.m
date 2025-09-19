function [Hc, Hc_normalized, Hcell] = compute_Hc_Hcell(BPs)
[~, nBase] = size(BPs);
Hcell = cell(1, nBase);
for iBase = 1:nBase
    Hcell{iBase} = ind2vec(BPs(:, iBase)')';
end

Hc = cell2mat(Hcell);
Hc_normalized = bsxfun(@rdivide, Hc, sqrt(sum(Hc, 1)));
end