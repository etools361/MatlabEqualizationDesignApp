function [corrected_WT2, corrected_AT2] = auto_repair_data(WT2, AT2)
    % 初始化参数
    [rows, cols] = size(WT2);
    corrected_WT2 = WT2;
    corrected_AT2 = AT2;
    
    % ================== 第一部分：检测异常列 ==================
    % 计算WT2每列残差并标记异常列
    wt2_residuals = zeros(size(WT2));
    p_wt2 = zeros(cols,3);
    for col = 1:cols
        x = (1:rows)';
        y = WT2(:, col);
        p_wt2(col,:) = polyfit(x, y, 2);
        wt2_residuals(:,col) = abs(y - polyval(p_wt2(col,:),x));
    end
    wt2_std = std(wt2_residuals);
    wt2_std(ismember(wt2_std,maxk(wt2_std,2))) = []; % 去除方差最大的两列
    wt2_threshold = max(wt2_std)*3.3;
    i_wt2 = find(std(wt2_residuals) > wt2_threshold);

    % 计算AT2每列残差并标记异常列
    at2_residuals = zeros(size(AT2));
    p_at2 = zeros(cols,2);
    for col = 1:cols
        x = (1:rows)';
        y = AT2(:, col);
        p_at2(col,:) = polyfit(x, y, 1);
        at2_residuals(:,col) = abs(y - polyval(p_at2(col,:),x));
    end
    at2_std = std(at2_residuals);
    at2_std(ismember(at2_std,maxk(at2_std,2))) = [];
    at2_threshold = max(at2_std)*3.3;
    i_at2 = find(std(at2_residuals) > at2_threshold);
    
    % ========== 第二部分：预测问题列的拟合系数 ==========
    % 预测WT2问题列系数（二次多项式）
    if ~isempty(i_wt2)
        prd_wt2_coeff = zeros(length(i_wt2),3);
        for k = 1:length(i_wt2)
            col = i_wt2(k);
            valid_cols = setdiff(1:cols,i_wt2);
            x_fit = valid_cols;
            for coeff = 1:3
                y_fit = p_wt2(valid_cols,coeff);
                p = polyfit(x_fit,y_fit,2); % 用二次曲线拟合系数变化
                prd_wt2_coeff(k,coeff) = polyval(p,col);
            end
        end
    end
    
    % 预测AT2问题列系数（一次多项式）
    if ~isempty(i_at2)
        prd_at2_coeff = zeros(length(i_at2),2);
        for k = 1:length(i_at2)
            col = i_at2(k);
            valid_cols = setdiff(1:cols,i_at2);
            x_fit = valid_cols;
            for coeff = 1:2
                y_fit = p_at2(valid_cols,coeff);
                p = polyfit(x_fit,y_fit,1); % 用直线拟合系数变化
                prd_at2_coeff(k,coeff) = polyval(p,col);
            end
        end
    end
    
    % ========== 第三部分：定位并修正异常点 ==========
    % 处理WT2异常点
    if ~isempty(i_wt2)
        for k = 1:length(i_wt2)
            col = i_wt2(k);
            % 计算预测值并找最大残差点
            y_pred = polyval(prd_wt2_coeff(k,:),1:rows);
            [~,max_idx] = max(abs(WT2(:,col)-y_pred'));
            % 重新拟合修正
            x_data = (1:rows)';
            y_data = corrected_WT2(:,col);
            y_data(max_idx) = NaN;
            valid = ~isnan(y_data);
            p_new = polyfit(x_data(valid),y_data(valid),2);
            corrected_WT2(max_idx,col) = polyval(p_new,max_idx);
        end
    end
    
    % 处理AT2异常点
    if ~isempty(i_at2)
        for k = 1:length(i_at2)
            col = i_at2(k);
            % 计算预测值并找最大残差点
            y_pred = polyval(prd_at2_coeff(k,:),1:rows);
            [~,max_idx] = max(abs(AT2(:,col)-y_pred'));
            % 重新拟合修正
            x_data = (1:rows)';
            y_data = corrected_AT2(:,col);
            y_data(max_idx) = NaN;
            valid = ~isnan(y_data);
            p_new = polyfit(x_data(valid),y_data(valid),1);
            corrected_AT2(max_idx,col) = polyval(p_new,max_idx);
        end
    end
end
