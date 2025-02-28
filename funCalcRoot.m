%--------------------------------------------------------------------------
% Edited by bbl
% Function: funCalcRoot
% Description: Calculate roots of a function based on given data and
% threshold.
%
% INPUTS:
%   t1:         Time or x data
%   y1:         Function values
%   varargin:
%     - Threshold:   Threshold value for roots
%     - NthRoot:     Specify Nth root (optional, default: 0 for all roots)
%     - SW_Interp:   Interpolation switch (optional, default: 1 for linear
%                    interpolation)
%
% OUTPUTS:
%   t:      Time or x values corresponding to the roots
%   rf:     Root flags indicating whether each root is a rising or falling
%           zero crossing (-1 for falling, 1 for rising)
%
% EXAMPLE:
%   [t, rf] = funCalcRoot(t1, y1, Threshold, NthRoot, SW_Interp);
%   [t, rf] = funCalcRoot(y1, Threshold);
%--------------------------------------------------------------------------
function [t, rf] = funCalcRoot(t1, y1, varargin)
    if nargin < 3
        t1T = t1;
        y1T = y1;
        t1 = 1:length(t1T);
        y1 = t1T;
        Threshold = y1T;
    else
        Threshold = varargin{1};
    end
    
    if nargin > 3
        NthRoot = varargin{2};
    else
        NthRoot = 0;
    end
    
    if nargin > 4
        SW_Interp = varargin{3};
    else
        SW_Interp = 1;
    end
    
    t = [];
    rf = [];
    
    [ay1, by1] = size(y1);
    
    % Check data format
    if (ay1 == 1) && (by1 >= 2)
    elseif (ay1 >= 2) && (by1 == 1)
        y1 = y1';
        t1 = t1';
    else
        warning(sprintf('Data y1 format error!\n'));
    end
    
    y2 = y1 - Threshold;
    sy2 = sign(y2);
    sy2(sy2 == 0) = 1; % Ensure sign change at zero crossings
    dsy2 = [diff(sy2), 0];
    [idsy2] = find(dsy2 ~= 0);
    [ifall] = find(dsy2 == -2);
    mRoots = length(idsy2);
    rf = ones(1, mRoots);
    rf(ismember(idsy2, ifall)) = -1;
    aidsy2 = [idsy2; idsy2 + 1];
    
    if mRoots
        if NthRoot
            if mRoots < NthRoot
                t = [];
                fprintf('Can''t find the Roots!\n');
                fprintf('mRoots=%d < NthRoot=%d\n', mRoots, NthRoot);
            else
                if SW_Interp
                    % Interpolate
                    t = interp1(y1(aidsy2(:, NthRoot)), t1(aidsy2(:, NthRoot)), Threshold);
                else
                    t = t1(idsy2(NthRoot));
                end
            end
        else
            % NthRoot = 0; for all roots.
            if SW_Interp
                % Interpolate for all roots
                for i = 1:mRoots
                    t(i) = interp1(y1(aidsy2(:, i)), t1(aidsy2(:, i)), Threshold);
                end
            else
                t = t1(idsy2);
            end
        end
    else
        t = [];
        fprintf('Can''t find the Roots!\n');
    end
end
