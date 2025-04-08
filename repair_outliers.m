function [WT_repaired, AT_repaired] = repair_outliers(WT, AT)
    % 修复WT矩阵
    WT_repaired = zeros(size(WT));
    for i = 1:size(WT, 1)
        y = WT(i, :);
        x = 1:length(y);
        [repaired_y, ~] = ransac_polyfit(x, y);
        WT_repaired(i, :) = repaired_y;
    end
    
    % 修复AT矩阵
    AT_repaired = zeros(size(AT));
    for i = 1:size(AT, 1)
        y = AT(i, :);
        x = 1:length(y);
        [repaired_y, ~] = ransac_polyfit(x, y);
        AT_repaired(i, :) = repaired_y;
    end
end

function [repaired_y, best_order] = ransac_polyfit(x, y)
    best_order = 1;
    repaired_y = y;
    min_total_residual = Inf;
    max_order = 3;
    max_iterations = 100;
    num_points = length(x);
    
    for order = 1:max_order
        required_samples = order + 1;
        if num_points < required_samples
            continue;
        end
        
        best_inliers = [];
        best_model = [];
        
        for iter = 1:max_iterations
            samples = randperm(num_points, required_samples);
            x_sample = x(samples);
            y_sample = y(samples);
            
            try
                p = polyfit(x_sample, y_sample, order);
            catch
                continue;
            end
            
            y_fit = polyval(p, x);
            residuals = abs(y - y_fit);
            mad_residual = mad(residuals, 1); % 中位数绝对偏差
            if mad_residual == 0
                mad_residual = eps;
            end
            inliers = find(residuals < 3 * mad_residual);
            
            if length(inliers) > length(best_inliers)
                best_inliers = inliers;
                best_model = p;
            end
        end
        
        if ~isempty(best_inliers) && length(best_inliers) >= required_samples
            x_in = x(best_inliers);
            y_in = y(best_inliers);
            p_final = polyfit(x_in, y_in, order);
            y_fit_all = polyval(p_final, x);
            residuals_all = abs(y - y_fit_all);
            sigma = std(residuals_all(best_inliers));
            if sigma == 0
                sigma = eps;
            end
            outliers = find(residuals_all >= 3 * sigma);
            
            y_repaired = y;
            y_repaired(outliers) = y_fit_all(outliers);
            valid_indices = setdiff(1:num_points, outliers);
            total_residual = sum(residuals_all(valid_indices));
            
            if total_residual < min_total_residual
                min_total_residual = total_residual;
                best_order = order;
                repaired_y = y_repaired;
            end
        end
    end
end