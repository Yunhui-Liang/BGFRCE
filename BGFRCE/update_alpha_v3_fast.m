function alpha = update_alpha_v3_fast(Gcell, HcZsqM, YinvYTY, Y)

nBase = length(Gcell);

% ==== 预计算 ====
P = YinvYTY * Y';   % n×n
Qcell = cell(nBase,1);
for j = 1:nBase
    Qcell{j} = HcZsqM' * Gcell{j};   % d×n
end

E = zeros(nBase);

% ==== 主循环 ====
parfor i = 1:nBase
    Ei = zeros(1,nBase);
    for j = i:nBase
        tmp1 = Qcell{j};          % d×n
        tmp2 = tmp1 - tmp1 * P;   % d×n
        tmp3 = tmp2 * Gcell{i} * HcZsqM;  % d×d
        Ei(j) = trace(tmp3);
    end
    E(i,:) = Ei;
end

% ==== 对称填充 ====
E = E + triu(E,1)';   % 保持对称

% ==== QP 求解 ====
H = (E + E');     % double-check symmetric
f = zeros(nBase,1);
lb = zeros(nBase,1);
ub = [];
Aeq = ones(1,nBase);
beq = 1;
opt = optimoptions('quadprog','Display','off');

alpha = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opt);

end
