function p = my_poly(roots)
    % 输入：roots 是根的数组，例如 [r1, r2, r3]
    % 输出：p 是多项式系数数组，例如 [1, a1, a2, ...]

    p = 1; % 初始化为 1（多项式 1）
    for i = 1:length(roots)
        p = my_conv(p, [1, -roots(i)]); % 逐步乘以 (s - r_i)
    end
end