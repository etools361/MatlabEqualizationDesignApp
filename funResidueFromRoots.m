%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 由零极点计算多项式的有理展开
%--------------------------------------------------------------------------
function [residues, poles, constant] = funResidueFromRoots(zeros, poles, num_lead, den_lead)
    % 输入参数增加分子分母最高次项系数（默认首一多项式）
    if nargin < 4
        den_lead = 1;
    end
    if nargin < 3
        num_lead = 1;
    end
    
    zeros = zeros(:).';
    poles = poles(:).';
    assert(length(zeros) == length(poles), '零极点数量必须相等');
    
    n = length(poles);
    residues = zeros(1, n);
    K = num_lead / den_lead; % 最高次项系数比
    
    for i = 1:n
        current_pole = poles(i);
        numerator = prod(current_pole - zeros);
        denominator = prod(current_pole - poles([1:i-1, i+1:n]));
        residues(i) = K * numerator / denominator;
    end
    
    constant = K;
end