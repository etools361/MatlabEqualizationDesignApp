N = 3;      % 项数
% m = 0.025;      % 区间左端点
% n = 0.55;      % 区间右端点
c = 10;    % 目标中心值
tol = 1e-6; % 收敛容差
max_iter = 50; % 最大迭代次数

% 调用函数
% [delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter);
kT = linspace(1.1,10,100);
nn = length(kT);
deltaT=[];
AT    =[];
WT    =[];
for ii=1:nn
    A_set = A;
    W_set = W;
    delta_set = delta;
    k = kT(ii);
    p = 1/k*0.9;
    m = 1/k*p;      % 区间左端点
    n = k*p;      % 区间右端点
    [delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, A_set, W_set,delta_set, 0);
    deltaT(ii) = delta;
    AT(ii,:) = A;
    WT(ii,:) = W;
    plot(kT(1:ii), deltaT(1:ii), '-r', 'linewidth', 2);
    grid on;
    drawnow;
end
AT1 = AT(:,1);
AT2 = AT(:,2);
AT3 = AT(:,3);


