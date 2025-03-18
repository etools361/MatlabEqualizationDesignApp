function y = my_polyval(p, x)
    % 输入：p 是多项式系数数组，x 是点的值
    % 输出：y 是多项式在 x 处的值

    y = 0;
    for i = 1:length(p)
        y = y + p(i) * x^(length(p) - i);
    end
end