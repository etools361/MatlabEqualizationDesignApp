%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 获取全部参数，由对称性获取参数
%--------------------------------------------------------------------------
% Helper function: Returns the current parameters a, w, and delta from pSol
function [aCurr, wCurr, deltaCurr] = funGetFullPara(pSol, W0, N)
    nn = floor((N+1)/2);
    if mod(N, 2)
        aCurr = [flip(pSol(2:nn)); pSol(1:nn)];
        wCurr = [flip(W0^2./pSol(nn+2:end-1)); pSol(nn+1:end-1)];
    else
        aCurr = [flip(pSol(1:nn)); pSol(1:nn)];
        wCurr = [flip(W0^2./pSol(nn+1:end-1)); pSol(nn+1:end-1)];
    end
    deltaCurr = pSol(end);
end