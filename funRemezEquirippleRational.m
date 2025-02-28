%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 使用Remez算法实现参数的计算
% REMEZ_EQUIRIPPLE_RATIONAL Uses a method based on the Remez exchange algorithm
% to compute the equiripple rational approximation of a rational function y_N(x)
% in the interval [m, n], aiming for the function to alternately equal c +/- delta 
% at 2N+1 points, with endpoints at c-delta.

% Inputs:
%   N       : The degree of the rational approximation (or order)
%   m, n    : The endpoints of the interval (m > 0, n > m)
%   c       : The middle value of the desired equiripple target
%   maxIter : Maximum number of iterations
%   tol     : Convergence threshold (for changes in delta or x_i)

% Outputs:
%   aSol, wSol, deltaSol : Converged parameters
%   xPointsSol          : Converged set of N+1 extrema points (including endpoints)

% Example usage:
%   >> [a,w,delta,xp] = remez_equiripple_rational(2, 1, 3, 5, 50, 1e-8);
%--------------------------------------------------------------------------
function [aSol, wSol, deltaSol, xPointsSol] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, dispEn, aSet, wSet)
    if nargin < 5
        maxIter = 20; % Default max iterations
    end
    if nargin < 6
        tol = 1e-8; % Default tolerance
    end
    if nargin < 7
        dispEn = 0; 
    end

    % -------------------------------
    % 1) Generate initial guess for N+1 points using logarithmic spacing
    % -------------------------------
    xLog = logspace(log10(sqrt(m*n)), log10(n), N+1);
    xPoints = xLog(:);  % Column vector of points
    W0 = sqrt(m*n);

    % -------------------------------
    % 2) Set initial guesses for (a_i, w_i, delta)
    % -------------------------------
    nn = floor((N+1)/2); % Half the number of points
    p0 = zeros(2*nn+1, 1); % Initial parameter vector    
    wh = n;wl = m;
    if N<3
        Ac = c;
    else
        K = c*log10(wh/wl)/(20*(N-1));
        Ac = 20*(10^K-1)/(10^K+1);
    end
    if nargin < 8
        p0(1:nn) = Ac;  % Initial guess for a
    else
        p0(1:nn) = aSet;  % Initial guess for a
    end
%     p0(1:nn) = aSet;  % Initial guess for a
%     pall = logspace(log10(m), log10(n), 2*N+1);
    R = log10(wh/wl);
    fRN = sqrt(wh*wl);
    if N~=1
        fRN0 = R/N*(0.4444*N^2 + 4.4*N -4.844) / (N + 9.364)+(1.587*N -1.213) / (N + 31.74);
%         fRN0 = R/N*(124.1*N-123.2)/(N+246.4)+(4.644*N)/(N+134);
        fRN0 = 10^fRN0;
    else
        fRN0 = 1;
    end
%     N
    fRN = fRN*fRN0;
%     Ac
    if mod(N,2)
        pall = logspace(log10(sqrt(wh*wl)), log10(fRN), floor((N+1)/2));
    else
        pall = logspace(log10(sqrt(1/fRN)), log10(fRN), floor((N+1)/2)*2);
    end
    kk = sqrt((1+Ac/20)/(1-Ac/20));
    mm = kk^(40/c).*1.09;%((1+Ac/20)/(1-Ac/20))^(20/c)
    if nargin < 9
%         if mod(N,2)
%             p0(nn+1:2*nn) =  W0.*mm.^(0:nn-1);  % Initial guess for w
%         else
%             p0(nn+1:2*nn) =  W0.*mm.^((0:nn-1)+0.5);  % Initial guess for w
%         end
        p0x = fliplr(pall);
        p0(nn+1:2*nn) = fliplr(p0x(1:nn))';  % Initial guess for w
    else
        p0(nn+1:2*nn) = wSet;
        % GEN x
        [A, W, delta] = funGetFullPara(p0, W0, N);
        xi = logspace(log10(sqrt(m*n)), log10(n), 10000);
        yVal = funCalcY(A, W, xi);
        dy = diff(yVal);
        [t, rf] = funCalcRoot(xi, dy, 0);
        xPoints0 = [sqrt(m*n),t,n];
        if length(xPoints0) == N+1 && N>5
            xPoints = xPoints0;
        end
    end
%     p0(nn+1:2*nn) = wSet;  % Initial guess for w
	rip = funGetRipple(Ac, mm);
    p0(end) = rip;    % Initial guess for delta (ensure >0)
    
    % -------------------------------
    % 3) Iteration: Exchange extrema points <-> Solve system
    % -------------------------------
    oldDelta = Inf; % Initialize previous delta to a large value
    for iter = 1:maxIter
        % (a) Solve for (a, w, delta) such that y_N(x_i) = c +/- delta
        [pSol, resNorm, exitflag] = solveSystemRemez(p0, xPoints, c, N, W0);
        if exitflag <= 0
            warning('fsolve may not have converged in this iteration.');
        end
        [aCurr, wCurr, deltaCurr] = funGetFullPara(pSol, W0, N);
        if bitand(dispEn,2)
            % Plot the result for visualization
            x = logspace(log10(xPoints(1)^2/xPoints(end)), log10(xPoints(end)), 1000);
            y = funCalcY(aCurr, wCurr, x);
            xPointsNew2 = findNewExtrema(aCurr, wCurr, c, m, n, N+1);
            semilogx(x, y, '-r');
            grid on;
            hold on;
            iy = interp1(x, y, xPointsNew2);
            semilogx(xPointsNew2, iy, '*r');
            hold off;
            title(sprintf('N:%d,Iter:%d', N,iter));
            xlabel('w/rand/s');
            ylabel('dH/dw');
            drawnow;
            pause(1);
        end
        % (b) Check for convergence based on delta change
        if abs(deltaCurr - oldDelta) < tol
            if bitand(dispEn,4)
                fprintf('Converged for delta in iteration %d: %.6e\n', iter, deltaCurr);
            end
            break;
        end

        % (c) Refine extrema points in [m, n]
        xPointsNew = findNewExtrema(aCurr, wCurr, c, m, n, N+1);

        % (d) Check for convergence based on extrema points change
        diffX = max(abs(xPointsNew - xPoints));
        if diffX < tol
            if bitand(dispEn,4)
                fprintf('Converged for x_i in iteration %d, max change = %.3e\n', iter, diffX);
            end
            xPoints = xPointsNew;
            break;
        end

        % (e) Update parameters for the next iteration
        xPoints = xPointsNew;
        p0 = pSol;
        oldDelta = deltaCurr;
        if bitand(dispEn,4)
            fprintf('Iteration %d complete: delta=%.6e, resNorm=%.3e, max dX=%.3e\n', ...
                     iter, deltaCurr, resNorm, diffX);
        end
    end

    % -------------------------------
    % Final output
    % -------------------------------
    aSol = pSol(1:nn);
    wSol = pSol(nn+1:2*nn);
    deltaSol = pSol(end);
    xPointsSol = xPoints;
    if bitand(dispEn,1)
        fprintf('=== Remez Equiripple Method Convergence Results ===\n');
        for i = 1:nn
            fprintf('  A_%d = %.6g\n', i, aSol(i));
        end
        for i = 1:nn
            fprintf('  W_%d = %.6g\n', i, wSol(i));
        end
        fprintf('  delta = %.6g\n', deltaSol);
%         disp('  x_i = ');
%         disp(xPointsSol(:).');
    end

end

% Helper function: Solves the system of equations using fsolve
function [pSol, resNorm, exitflag] = solveSystemRemez(p0, xPoints, c, N, W0)
    eqFun = @(p) systemEquations(p, xPoints, c, N, W0);
    options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
    [pSol, ~, exitflag, output] = fsolve(eqFun, p0, options);
    Fval = eqFun(pSol);
    resNorm = norm(Fval, 2);
end

% Helper function: Constructs the system of equations for fsolve
function F = systemEquations(p, xPoints, c, N, W0)
    if mod(N, 2)
        nn = floor((N+1)/2);
        p(nn+1) = W0;
    end
    [A, W, delta] = funGetFullPara(p, W0, N);
    nPts = length(xPoints); % Should be N+1
    F = zeros(nPts, 1);

    for i = 1:nPts
        xi = xPoints(i);
        yVal = funCalcY(A, W, xi);
        if mod(N, 2)
            if mod(i, 2) == 0
                F(i) = yVal - (c - delta);
            else
                F(i) = yVal - (c + delta);
            end
        else
            if mod(i, 2) == 0
                F(i) = yVal - (c + delta);
            else
                F(i) = yVal - (c - delta);
            end
        end
    end
end

% Helper function: Finds new extrema points within the interval [m, n]
function xPointsNew = findNewExtrema(A, W, c, m, n, numPoints)
    nFine = 500*length(numPoints);
    xFine = logspace(log10(sqrt(m*n)), log10(n), nFine);
    errFine = funCalcY(A, W, xFine) - c;
    errFineDeriv = diff(errFine) ./ diff(xFine);

    xExtrema = [sqrt(m*n); n];
    for j = 2:(nFine-1)
        if (errFineDeriv(j-1) > 0 && errFineDeriv(j) < 0) || ...
           (errFineDeriv(j-1) < 0 && errFineDeriv(j) > 0)
            xExtrema = [xExtrema; xFine(j)];
        end
    end

    xExtrema = unique(sort(xExtrema));
    shortNum = numPoints - length(xExtrema);
    if shortNum > 0
        xExtrema = fillLargestLogGap(xExtrema, shortNum);
    end

    if length(xExtrema) > numPoints
        if mod(numPoints-1, 2)
            logx = log10(xExtrema);
            M = diff(logx);
            [~, iM] = min(M);
            if iM == 1
                xExtrema(iM+1) = [];
            end
        else
            while length(xExtrema) > numPoints
                logx = log10(xExtrema);
                M = diff(logx);
                [~, iM] = min(M);
                xExtrema(iM) = [];
            end
        end
    end

    xPointsNew = sort(xExtrema);
end

% Helper function: Fills the largest log gap with new points
function xPointsNew = fillLargestLogGap(xExtrema, shortNum)
    if shortNum <= 0
        xPointsNew = xExtrema;
        return;
    end
    xExtrema = sort(unique(xExtrema));
    nExt = length(xExtrema);
    if nExt < 2
        xmin = min(xExtrema);
        xmax = max(xExtrema);
        if nExt == 1
            xmax = 2 * xmin;
        end
        newPoints = logspace(log10(xmin), log10(xmax), shortNum+2).';
        newPoints([1, end]) = [];
        xPointsNew = unique(sort([xExtrema; newPoints]));
        return;
    end

    logX = log(xExtrema);
    dlog = diff(logX);
    [maxVal, idxMax] = max(dlog);

    leftPt = xExtrema(idxMax);
    rightPt = xExtrema(idxMax+1);

    newPoints = logspace(log10(leftPt), log10(rightPt), shortNum+2).';
    newPoints([1, end]) = [];

    xPointsNew = unique(sort([xExtrema; newPoints]));
end


