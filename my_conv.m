function R = my_conv(P, Q)
    % 输入：P 和 Q 是两个多项式的系数数组
    % 输出：R 是卷积结果的多项式系数数组

    % 初始化结果数组
    len_P = length(P);
    len_Q = length(Q);
    R = zeros(1, len_P + len_Q - 1);

    % 计算卷积
    for i = 1:len_P
        for j = 1:len_Q
            R(i + j - 1) = R(i + j - 1) + P(i) * Q(j);
        end
    end
end