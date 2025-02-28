%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 三次样条插值计算
%--------------------------------------------------------------------------
function yq = funSplineInterp(x, y, xq)
    % 输入x和y是行向量，xq是查询点数组，输出yq为插值结果
    x = x(:).';
    y = y(:).';
    n = length(x) - 1; % 区间数
    h = x(2:n+1) - x(1:n); % 步长数组

    if any(h <= 0)
        error('x必须严格递增');
    end

    % 构建右侧向量d
    m = n - 1;
    d = zeros(1, m);
    for k = 1:m
        term1 = (y(k+2) - y(k+1)) / h(k+1);
        term2 = (y(k+1) - y(k)) / h(k);
        d(k) = 6 * (term1 - term2);
    end
    d = d(:);

    % 构建三对角矩阵A
    main_diag = 2 * (h(1:m) + h(2:m+1));
    lower_diag = h(2:m);
    upper_diag = h(2:m);

    % 应用Thomas算法求解AM = d
    a = [0; lower_diag(:)];
    b = main_diag(:);
    c = [upper_diag(:); 0];
    d_vec = d;
    m_solve = m;

    c_prime = zeros(m_solve, 1);
    d_prime = zeros(m_solve, 1);
    c_prime(1) = c(1) / b(1);
    d_prime(1) = d_vec(1) / b(1);

    for i = 2:m_solve-1
        denominator = b(i) - a(i) * c_prime(i-1);
        c_prime(i) = c(i) / denominator;
        d_prime(i) = (d_vec(i) - a(i)*d_prime(i-1)) / denominator;
    end

    denominator = b(m_solve) - a(m_solve)*c_prime(m_solve-1);
    d_prime(m_solve) = (d_vec(m_solve) - a(m_solve)*d_prime(m_solve-1)) / denominator;

    M = zeros(m_solve, 1);
    M(m_solve) = d_prime(m_solve);
    for i = m_solve-1:-1:1
        M(i) = d_prime(i) - c_prime(i) * M(i+1);
    end

    % 构建完整的M向量
    M_full = [0; M; 0];

    % 计算每个区间的系数
    coeffs = zeros(n, 4);
    for k = 1:n
        a_k = y(k);
        h_k = h(k);
        M_k = M_full(k);
        M_k_plus_1 = M_full(k+1);
        c_k = M_k / 2;
        d_k = (M_k_plus_1 - M_k) / (6 * h_k);
        b_k = (y(k+1) - y(k)) / h_k - h_k * (2*M_k + M_k_plus_1) / 6;
        coeffs(k, :) = [a_k, b_k, c_k, d_k];
    end

    % 对每个xq进行插值
    yq = zeros(size(xq));
    for i = 1:numel(xq)
        xi = xq(i);
%         if xi < x(1) || xi > x(end)
%             error('查询点超出范围');
%         end
        k = find(x <= xi, 1, 'last');
        if isempty(k)
            k = 1;
        end
        if k == n+1
            k = n;
        end
        dx = xi - x(k);
        a = coeffs(k, 1);
        b = coeffs(k, 2);
        c = coeffs(k, 3);
        d = coeffs(k, 4);
        yq(i) = a + b*dx + c*dx^2 + d*dx^3;
    end
end