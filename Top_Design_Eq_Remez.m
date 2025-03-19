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


