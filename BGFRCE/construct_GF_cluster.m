function [M, M_sqrt] = construct_GF_cluster(Hc, k, order)
[nFea] = size(Hc, 2);
options = [];
options.NeighborMode = 'KNN';
options.k = k;
options.WeightMode = 'HeatKernel';
S = constructW(Hc', options);
% M = compute_graph_power_mul(S, order);
% 适合高阶
d = max(sum(S, 2), eps);
DS = bsxfun(@times, S, d.^-0.5);
DSD = bsxfun(@times, DS, d'.^-0.5);
S = 0.5*(speye(nFea) + DSD);
S = (S + S')/2;
[M, M_sqrt] = compute_graph_power_eig_v3(S, order);

end