function [label_Y, objHistory] = BGFRCE(Hc, Hcell, label_Y, nCluster, k, order)

% nBase = length(Hcell);
[nSmp, ~] = size(Hc);
nBase = length(Hcell);
alpha = ones(1, nBase)./nBase;

% construct Graph filter G
Gcell = construct_GF_cell(Hcell);

[~, M_sqrt] = construct_GF_cluster(Hc, k, order);

objHistory = [];
maxIter = 20;
% compute G_star
G_star = sparse(nSmp);
for i1 = 1:nBase
    G_star = G_star + alpha(i1) * Gcell{i1};
end
G_star = (G_star + G_star')/2;
for iter = 1:maxIter
    %**************************
    % update Y
    %**************************
    HcsqM = Hc * M_sqrt;
    GHcsqM = G_star * HcsqM;
    
    
    [label_Y, ~, ~] = CDKM_fast(GHcsqM', label_Y, nCluster);
    
    % compute R
    Y = ind2vec(label_Y')';
    % Use pseudo-inverse for numerical stability
    YinvYTY = Y./sum(Y, 1);
    R = speye(nSmp) - YinvYTY * Y';% speye for sparse identity
    
    
    %**************************
    % update alpha
    %**************************
    %     alpha = update_alpha(Gcell, HcZsqM, R);
    
    alpha = update_alpha_v3_fast(Gcell, HcsqM, YinvYTY, Y);
    % compute G_star
    G_star = sparse(nSmp);
    for i1 = 1:nBase
        G_star = G_star + alpha(i1) * Gcell{i1};
    end
    
    obj = compute_obj(G_star, Hc, M_sqrt, R);
    objHistory = [objHistory; obj]; %#ok
    
    if iter > 2 && abs((objHistory(end-1) - objHistory(end)) / objHistory(end-1)) < 1e-6
        break;
    end
    
end
end



function obj = compute_obj(G_star, Hc, M_sqrt, R)

G_starHcM_sqrt = G_star * Hc * M_sqrt;
KR = G_starHcM_sqrt' * R * G_starHcM_sqrt;
KR = (KR + KR')/2;
obj = trace(KR);
end
