function alpha = update_alpha(Gcell, HcZsqM, R)


% HcZMZHc = HcZsqM * HcZsqM';
nBase = length(Gcell);
Q = zeros(nBase);
for i = 1:nBase
    for j = 1:nBase
        tmp1 = Gcell{i}*HcZsqM;
        tmp2 = HcZsqM' * Gcell{j};
        tmp3 = tmp1*tmp2 * R;
        Q(i,j) = sum(sum(tmp3.^2));
    end
end

H = (Q + Q');
f = zeros(nBase,1);

lb = zeros(nBase,1);
ub = [];

Aeq = ones(1,nBase);      % alpha sum = 1
beq = 1;

% opt = optimset('quadprog');
% opt.Algorithm = 'interior-point-convex';
% opt.Display = 'off';

opt = optimoptions('quadprog','Display','off');

[alpha,~, ~, ~, ~] = quadprog(H, f, [], [], Aeq, beq, lb, ub,[],opt);
end