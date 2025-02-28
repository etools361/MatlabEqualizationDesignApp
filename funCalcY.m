%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 计算斜率函数值
%--------------------------------------------------------------------------
% Helper function: Computes y_N(x)
function val = funCalcY(A, W, x)
    xx = x.^2;
    val = 0;
    k = sqrt((1 + A ./ 20) ./ (1 - A ./ 20));
    a = W .* k;
    b = W ./ k;
    for i = 1:length(a)
        val = val + 20 .* xx .* (a(i).^2 - b(i).^2) ./ ((xx + a(i).^2) .* (xx + b(i).^2));
    end
end