function Y0 = funinterp1(Xin, Yin, X0)
    % Find maximum and minimum x-values and their indices
    [XinMax, XinMaxIndex] = max(Xin);
    [XinMin, XinMinIndex] = min(Xin);
    mXin = length(Xin);
    mX0 = length(X0);
    for ii=1:mX0
        % Check if X0 is inside or outside the range
        if (XinMax - X0(ii)) * (X0(ii) - XinMin) >= 0 % Inside
            % Handle special cases when length of Xin is less than 3
            if mXin < 3
                Y0(ii) = GetY0(Xin, Yin, X0(ii));
            else
                % Check for repeated x-values and remove duplicates
                x_diff0 = diff(Xin);
                [ReX, ReXIndex] = find(x_diff0 == 0);
                if ~isempty(ReX)
                    Xin(ReXIndex) = [];
                    Yin(ReXIndex) = [];
                end

                x_diff = sum(diff(x_diff0));
                Xin1 = Xin(1:end - 1);
                Xin2 = Xin(2:end);
                x_div = sum(diff(Xin1 ./ Xin2));

                % Check if differences are zero and perform interpolation
                if or(x_diff == 0, x_div == 0)
                    Y0(ii) = interp1(Xin, Yin, X0(ii));
                else
                    % Handle edge cases where X0 is equal to an x-value
                    XX1 = Xin - X0(ii);
                    XX2 = sign(XX1);
                    IndexEq = find(XX2 == 0, 1);
                    if isempty(IndexEq)
                        XX3 = diff(XX2);
                        IndexCG = find(XX3 ~= 0, 1);
                        Y0(ii) = GetY0([Xin(IndexCG:IndexCG + 1)], Yin(IndexCG:IndexCG + 1), X0(ii));
                    else
                        Y0(ii) = Yin(IndexEq);
                    end
                end
            end
        else % Outside
            Xin1 = Xin;
            Yin1 = Yin;
            if X0(ii) > XinMax
                % Handle cases where X0 is greater than maximum x-value
                if XinMaxIndex > XinMinIndex
                    if XinMaxIndex > 2
                        Xin1 = Xin(XinMaxIndex - 1:XinMaxIndex);
                        Yin1 = Yin(XinMaxIndex - 1:XinMaxIndex);
                    end
                else
                    if mXin >= XinMaxIndex + 1
                        Xin1 = Xin(XinMaxIndex:XinMaxIndex + 1);
                        Yin1 = Yin(XinMaxIndex:XinMaxIndex + 1);
                    end
                end
            else % X0 is less than minimum x-value
                if XinMaxIndex > XinMinIndex
                    if mXin >= XinMinIndex + 1
                        Xin1 = Xin(XinMinIndex:XinMinIndex + 1);
                        Yin1 = Yin(XinMinIndex:XinMinIndex + 1);
                    end
                else
                    if XinMinIndex > 2
                        Xin1 = Xin(XinMinIndex - 1:XinMinIndex);
                        Yin1 = Yin(XinMinIndex - 1:XinMinIndex);
                    end
                end
            end
            % Perform interpolation for points outside the range
            Y0(ii) = GetY0(Xin1, Yin1, X0(ii));
        end
    end
end

% Helper function to calculate Y0 based on linear interpolation
function Y0 = GetY0(Xin, Yin, X0)
    if Xin(1) ~= Xin(2)
        k = (Yin(1) - Yin(2)) / (Xin(1) - Xin(2));
        b = Yin(1) - k .* Xin(1);
        Y0 = k .* X0 + b;
    else
        Y0 = [];
        warning('x1(%f) == x2(%f)\n', Xin(1), Xin(2));
    end
end
