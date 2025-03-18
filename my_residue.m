function [coeffs, poles, k] = my_residue(P, poles)
    % 输入：
    % P: 分子多项式系数数组，例如 [1, 2, 3] 表示 s^2 + 2s + 3
    % poles: 已知的极点数组，例如 [p1, p2, p3]
    % 输出：
    % coeffs: 部分分式系数
    % poles: 极点（直接返回输入）
    % k: 直接项（如果分子的次数大于或等于分母的次数）

    % 初始化部分分式系数
    coeffs = zeros(size(poles));

    % 构造分母多项式 Q
    Q = my_poly(poles);

    % 计算直接项（如果分子的次数大于或等于分母的次数）
    n = length(P) - 1; % 分子的次数
    m = length(Q) - 1; % 分母的次数
    if n >= m
        % 多项式除法，计算直接项
        [k, P] = my_deconv(P, Q);
    else
        k = []; % 没有直接项
    end

    % 计算部分分式系数
    for i = 1:length(poles)
        % 计算分母的导数在极点处的值
        denom = 1;
        for j = 1:length(poles)
            if j ~= i
                denom = denom * (poles(i) - poles(j));
            end
        end
        % 计算部分分式系数
        coeffs(i) = my_polyval(P, poles(i)) / denom;
    end
end