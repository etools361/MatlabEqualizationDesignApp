%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 计算多项式的有理分解
%--------------------------------------------------------------------------
function residues = funComputeResiduesForRL(a, b, Rs)
% 确保a和b是行向量
a = a(:).';
b = b(:).';

N = length(a);
residues = zeros(1, N);

for i = 1:N
    % 计算p(b_i) = product_{k=1}^N (b_i - a_k)
    p_val = 1;
    for k = 1:N
        p_val = p_val * (b(i) - a(k));
    end
    
    % 计算q’(b_i) = product_{k≠i} (b_i - b_k)
    q_prime_val = 1;
    for k = 1:N
        if k ~= i
            q_prime_val = q_prime_val * (b(i) - b(k));
        end
    end
    
    % 计算留数
    residues(i) = (2 * p_val) / (Rs * q_prime_val);
end
end