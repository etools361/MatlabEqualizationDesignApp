function [a1, b1] = solve_a1_b1(A, w)
    % 定义方程组
    fun = @(uv) [
        3*w^4 + w^2*(uv(1) + uv(2)) - uv(1)*uv(2);
        20*w*(uv(1) - uv(2)) / ((w^2 + uv(1))*(w^2 + uv(2))) - A
    ];
    
    % 初始猜测（需根据实际情况调整）
    uv0 = [2*w^2; w^2]; % 示例：假设 a1 ≈ sqrt(2)w, b1 ≈ w
    
    % 使用fsolve求解
    options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt');
    uv = fsolve(fun, uv0, options);
    
    a1 = sqrt(uv(1));
    b1 = sqrt(uv(2));
end