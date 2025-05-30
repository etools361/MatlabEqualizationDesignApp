%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 使用Remez算法计算参数
%   N       : The degree of the rational approximation (or order)
%   m, n    : The endpoints of the interval (m > 0, n > m)
%   c       : The middle value of the desired equiripple target
%   maxIter : Maximum number of iterations
%   tol     : Convergence threshold (for changes in delta or x_i)
%--------------------------------------------------------------------------
%     
N = 4;
m = 0.01*2*pi;
n = 100*2*pi;
delta = 3;
c = delta/log10(n/m);
maxIter = 20;
tol = 1e-6;
[aSol,wSol,deltaSol,xSol] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, 7);
[A, W, delta] = funGetFullPara([aSol;wSol;0], sqrt(m*n), N);
Slope = 0;% 0:Positive Slope;1:Negative Slope
Type = 1;% netlist:0,zobel network, 1,RC/RL Serial, 2, RC/RL Parallel
funGenSchAndSim(A, W, Rs, Slope, Type, m/2/pi, n/2/pi);



% 输入参数
N = 7;      % 项数
m = 0.025;      % 区间左端点
n = 0.55;      % 区间右端点
c = 5;    % 目标中心值
tol = 1e-6; % 收敛容差
max_iter = 500; % 最大迭代次数

% 调用函数
[delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter);
% 0.01,0.026,0.056,0.084,0.118,0.21,0.264,0.797,0.945
% 2.40,1.653,1.533,1.824,1.564,2.04,1.847,4.389,4.329

A_set = [];
W_set = [];
delta_set = 0.001;
[delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, A_set, W_set,delta_set);

semilogx(W,A,'-*', 'linewidth', 2);grid on;
x = linspace(m,n,1000);
y = funCalcY_Linear((A),(W),x);
semilogx(x, y, '-r', 'linewidth', 2);
grid on;
hold on;
for ii=1:N
    y = funCalcY_Linear((A(ii)),(W(ii)),x);
    plot(x, y, '-', 'linewidth', 2);
%     yf = A(ii)-(A(ii).*(300-A(ii)^2*W(ii)^2)./(400.*W(ii)^3)).*(x-W(ii)).^2;
%     plot(x, yf, '-', 'linewidth', 2);
end
hold off;

% 显示结果
disp('Delta:'); disp(delta);
disp('A:'); disp(A);
disp('W:'); disp(W);
disp('极值点:'); disp(x_points);
A_set = A;
W_set = W;
delta_set = delta;
k = 1.414;
p = 1/k*0.9;
m = 1/k*p;      % 区间左端点
n = k*p;      % 区间右端点
[delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, A_set, W_set,delta_set);
[val, x] =  funEvalLinearEq(A, W, [], 0);
hold on;
plot(x, val, '-b', 'linewidth', 2);
hold off



