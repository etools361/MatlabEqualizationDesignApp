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
N = 3;
m = 0.01;
n = 100;
c = 10;
maxIter = 20;
tol = 1e-6;
[aSol,wSol,deltaSol,xSol] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, 7);



