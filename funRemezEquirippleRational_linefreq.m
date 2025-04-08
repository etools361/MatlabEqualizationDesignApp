function [delta, A, W, x_points, Err] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, a_guess, w_guess, delta_guess, dispEn)
    if nargin<7
        a_guess  = c.*(1:N)./N; % A的初始猜测
    end
    if nargin<8
        w_guess  = logspace(log10(m)*0.9, log10(n), N); % a的初始猜测
    end
    if nargin<9
        delta_guess = 0.1 * c/N; % 初始delta猜测
    end
    if nargin<10
        dispEn = 1;
    end
    xFix = 0;
    % 参数初始化
    x_points = logspace(log10(m), log10(n), 2*N+1);%1,3;2,5;3,7
    for iter = 1:max_iter
        p0 = [a_guess,w_guess,delta_guess];
        % 求解方程组以获得新的参数
        [pSol, resNorm, exitflag] = solveSystemRemez(p0, x_points, c, N);
        if exitflag <= 0
            warning('fsolve may not have converged in this iteration.');
        end    
        A = pSol(1:N);
        W = pSol(N+1:2*N);
        delta = pSol(2*N+1);
        x_samples = logspace(log10(m), log10(n), 1000);
        val = funCalcY_Linear(A, W, x_samples);
        if dispEn
            plot(x_samples, val, '-r', 'linewidth', 2);grid on;
    %         hold on;
    %         x_samples = linspace(0, 1, 1000);
    %         val = funCalcY_Linear(A, W, x_samples);
    %         plot(x_samples, val, '--b', 'linewidth', 2);grid on;
    %         hold off;
            drawnow;
        end
%         pause(1);
        
        % 寻找新的极值点
        [new_x, randx] = find_new_extremes(m, n, A, W, c, delta, N);
        if length(new_x)>2*N+1
            new_x = new_x(1:2*N+1);
        end
        if max(abs(new_x - x_points)) < tol
            if min(W.*A)<0
                fprintf('Randn x\n');
                if xFix == 0
                    xFix = 1;
                    W0 = abs(W);
                    new_x0 = new_x;
                else
                    W = abs(W0);
                    new_x = new_x0;
                end
                A     = c.*(1:N)./N;
            else
                [W, iW] = sort(abs(W));
                A = abs(A(iW));
                break;
            end
        end
%         if xFix == 0 
%             if randx
%                 W = logspace(log10(m)*0.9, log10(n), N);
%                 new_x = logspace(log10(m), log10(n), 2*N+1);
%                 A     = randi([10,20],1,1)./10.*c.*(1:N)./N;
%             end
%         else
%             if randx
%                 W = W0;
%                 new_x = new_x0;
%                 A     = randi([10,20],1,1)./10.*c.*(1:N)./N;
%             end
%         end
        
        % 更新极值点和猜测
        x_points = new_x;
        if delta<0
            delta = abs(delta);
        end
        delta_guess = delta;
        a_guess = A;
        w_guess = W;
    end
    Err = 0;
    if iter==max_iter
        Err = 1;
    end
end

function [pSol, resNorm, exitflag] = solveSystemRemez(p0, xPoints, c, N)
    eqFun = @(p) systemEquations(p, xPoints, c, N);
    options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
    [pSol, ~, exitflag, output] = fsolve(eqFun, p0, options);
%     [pSol, ~, exitflag] = my_newton_method(eqFun, p0, 1e-7, 200);
    %4.5167    2.8015    2.4127    2.3117    2.3074    2.3662    2.5423    3.0935    9.1160
    %0.0046    0.0106    0.0197    0.0348    0.0608    0.1068    0.1926    0.3766    1.9000
    Fval = eqFun(pSol);
    resNorm = norm(Fval, 2);
end

% Helper function: Constructs the system of equations for fsolve
function F = systemEquations(p, xPoints, c, N)
    A = p(1:N);
    W = p(N+1:2*N);
    delta = p(2*N+1);
    nPts = length(xPoints); % Should be N+1
    F = zeros(nPts, 1);

    for i = 1:nPts
        xi = xPoints(i);
        yVal = funCalcY_Linear(A, W, xi);
        if mod(i, 2) == 0
            F(i) = yVal - (c + delta);
        else
            F(i) = yVal - (c - delta);
        end
    end
end

function [new_x, randx] = find_new_extremes(m, n, a, b, c, delta, N)
    % 寻找导数为零的点
    x_samples = logspace(log10(m), log10(n), 1000);
    e = funCalcY_Linear(a, b, x_samples) - c;
    
    % 寻找极大值和极小值
    try
        [~, max_indices] = findpeaks(real(e));
        [~, min_indices] = findpeaks(-real(e));
    catch
        [~, max_indices] = findpeaks(e);
        [~, min_indices] = findpeaks(-e);
    end
    candidates = sort([x_samples(max_indices), x_samples(min_indices)]);
    % 加上端点
    candidates = [m,candidates,n];
    randx = 0;
    % 确保有足够的候选点
    NN = length(candidates);
    if NN < 2*N+1
        for ii=1:2*N+1-NN
            [maxCand, iCnad] = max(diff(log10(candidates)));
            candidates = [candidates(1:iCnad),sqrt(candidates(1+iCnad)*candidates(iCnad)),candidates(1+iCnad:end)];
        end
        randx = 1;
        new_x = candidates;
    else
        new_x = candidates;
    end
    
    % 按误差绝对值排序
%     e_candidates = funCalcY_Linear(a, b, candidates) - c;
%     [~, idx] = sort(abs(e_candidates), 'descend');
%     selected = candidates(idx(1:2*N+1));
    
    % 确保符号交替
%     e_sign = sign(e_candidates(idx(1:2*N+1)));
%     valid = true;
%     for i = 2:length(e_sign)
%         if e_sign(i) == e_sign(i-1)
%             valid = false;
%             break;
%         end
%     end
%     
%     if valid
%         new_x = sort([m; selected(2:end-1); n]); % 固定端点
%         success = true;
%     else
%         new_x = [];
%         success = false;
%     end
end

function [pSol, resNorm, exitflag] = my_newton_method(eqFun, p0, tol, max_iter)
    exitflag = 1; % 1-收敛, 0-奇异, -1-未收敛
    p = p0(:);    % 强制转为列向量
    n = length(p);
    lambda = 1e-6; % 初始正则化参数
    h_min = 1e-10; % 最小步长
    tol_ripple = 1e-6;
    
    for iter = 1:max_iter
        F = eqFun(p);
        resNorm = norm(F);
        current_ripple = max(F) - min(F); % 计算当前纹波
        
        % 双重收敛条件：残差或纹波达标
        if resNorm < tol || current_ripple < tol_ripple
            exitflag = 1; % 标记为收敛
            break;
        end
        
        % 动态步长：每个参数分量独立计算步长
        h = max(1e-8 * (1 + abs(p)), h_min); % 正确维度操作
        
        % 数值计算雅可比矩阵
        J = zeros(n, n);
        for i = 1:n
            p_perturbed = p;
            p_perturbed(i) = p_perturbed(i) + h(i); % 使用分量步长
            F_perturbed = eqFun(p_perturbed);
            J(:, i) = (F_perturbed - F) / h(i);
        end
        
        % 正则化处理
        J_reg = J' * J + lambda * eye(n);
        grad = J' * (-F);
        
        % 解方程：避免直接使用 \
        if rcond(J_reg) < eps
            exitflag = 0; % 标记为奇异
            warning('正则化矩阵仍接近奇异');
            break;
        end
        delta = J_reg \ grad;
        
        % 更新参数前检查残差
        p_new = p + delta;
        F_new = eqFun(p_new);
        resNorm_new = norm(F_new);
        
        if resNorm_new < resNorm
            p = p_new;
            lambda = max(lambda / 10, 1e-12); % 成功则减小正则化
        else
            lambda = lambda * 10; % 失败则增大正则化
            % 此处不更新 p，保持当前值
        end
    end
    
    % 处理未收敛情况
    if iter >= max_iter && resNorm >= tol
        exitflag = -1;
        warning('未在最大迭代次数内收敛');
    end
    
    pSol = p'; % 返回列向量
    resNorm = norm(eqFun(p));
end
