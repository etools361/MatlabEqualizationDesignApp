%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 计算斜率函数值
%--------------------------------------------------------------------------
function val = funCalcY_Linear(A, W, x)
    xx = x.^2;
    val = 0;
    aa = log(10);
    a = W .* sqrt((30/aa+A.*W)./(10/aa-A.*W));
    b = W .* sqrt((30/aa-A.*W)./(10/aa+A.*W));
    for i = 1:length(a)
        val = val + 20./aa .* x .* (a(i).^2 - b(i).^2) ./ ((xx + a(i).^2) .* (xx + b(i).^2));
    end
end
