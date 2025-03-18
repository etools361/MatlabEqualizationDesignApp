function [residues_r2, poles_r] = computeResidues(PP, QQ)
    % 确保a和b是行向量
%     a = a(:).';
%     b = b(:).';
% 
%     N = length(a);
%     residues_r2 = zeros(1, N);
%     poles_r = zeros(1, N);

    % 手动实现 poly 函数
%     Ps = myPoly(-b);
%     Qs = myPoly(-a);

    % 计算 P0 和 Q0
%     P0 = prod(b);
%     Q0 = prod(a);

    % 计算 PP 和 QQ
%     QQ = Qs;
%     PP = polySubtract(polyMultiply(Ps, Q0 / P0), Qs);
    PP0 = PP(1);

    % 计算留数和极点
    for i = 1:N
        % 计算 p(b_i) = product_{k=1}^N (b_i - a_k)
        p_val = 1;
        for k = 1:N
            p_val = p_val * (QQ(i) - a(k));
        end

        % 计算 q’(b_i) = product_{k≠i} (b_i - b_k)
        q_prime_val = 1;
        for k = 1:N
            if k ~= i
                q_prime_val = q_prime_val * (QQ(i) - QQ(k));
            end
        end

        % 计算留数
        residues_r2(i) = (p_val) / (q_prime_val * PP0);
        poles_r(i) = -QQ(i);
    end

end

function p = myPoly(roots)
    % 手动实现 poly 函数，根据根生成多项式系数
    N = length(roots);
    p = [1]; % 初始化多项式为 1

    for i = 1:N
        % 每次乘以 (s - roots(i))
        p = polyMultiply(p, [1, -roots(i)]);
    end
end

function result = polyMultiply(p, q)
    % 手动实现多项式乘法
    m = length(p) - 1; % p 的最高次数
    n = length(q) - 1; % q 的最高次数
    result = zeros(1, m + n + 1); % 初始化结果多项式

    for i = 0:m
        for j = 0:n
            result(i + j + 1) = result(i + j + 1) + p(i + 1) * q(j + 1);
        end
    end
end

function result = polySubtract(p, q)
    % 手动实现多项式减法
    m = length(p);
    n = length(q);
    maxLen = max(m, n);

    % 补齐长度
    p = [zeros(1, maxLen - m), p];
    q = [zeros(1, maxLen - n), q];

    result = p - q;
end