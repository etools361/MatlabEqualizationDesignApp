function [cleaned_WT3, cleaned_AT3] = jointRepairOutliers(WT3, AT3, polynomial_order, threshold)
    % 输入参数:
    % WT3, AT3 - 相关联的m*n二维数据数组
    % polynomial_order - 多项式阶数(1-3)
    % threshold - 联合残差阈值(建议3-4)
    
    % 初始化输出
    cleaned_WT3 = WT3;
    cleaned_AT3 = AT3;
    
    % 获取数据维度
    [m, n] = size(WT3);
    
    % 创建索引向量
    x = (1:m)';
    
    for col = 1:n
        % 当前列数据
        wt = WT3(:, col);
        at = AT3(:, col);
        
        % 方法1: 双变量联合多项式拟合
        % 拟合WT3趋势
        p_wt = polyfit(x, wt, polynomial_order);
        fit_wt = polyval(p_wt, x);
        
        % 拟合AT3趋势
        p_at = polyfit(x, at, polynomial_order);
        fit_at = polyval(p_at, x);
        
        % 计算标准化残差
        res_wt = (wt - fit_wt)/std(wt - fit_wt);
        res_at = (at - fit_at)/std(at - fit_at);
        
        % 方法2: 交叉验证回归（WT3和AT3互相预测）
        % 建立WT3->AT3回归模型
        p_cross = polyfit(wt, at, polynomial_order);
        pred_at = polyval(p_cross, wt);
        
        % 建立AT3->WT3回归模型
        p_cross_rev = polyfit(at, wt, polynomial_order);
        pred_wt = polyval(p_cross_rev, at);
        
        % 计算交叉残差
        cross_res_at = (at - pred_at)/std(at - pred_at);
        cross_res_wt = (wt - pred_wt)/std(wt - pred_wt);
        
        % 组合四种残差指标
        combined_res = sqrt(res_wt.^2 + res_at.^2 + cross_res_wt.^2 + cross_res_at.^2);
        
        % 检测异常点
        outliers = combined_res > threshold;
        
        % 三维曲面修复（考虑相邻数据点）
        for i = 1:m
            if outliers(i)
                % 使用滑动窗口局部拟合
                window_size = 7;
                window_start = max(1, i - window_size);
                window_end = min(m, i + window_size);
                valid_window = window_start:window_end;
                valid_window = valid_window(~outliers(valid_window)); % 排除窗口内异常点
                
                % 局部多项式拟合
                if length(valid_window) > 5
                    p_local_wt = polyfit(x(valid_window), wt(valid_window), polynomial_order);
                    repaired_wt = polyval(p_local_wt, x(i));
                    
                    p_local_at = polyfit(x(valid_window), at(valid_window), polynomial_order);
                    repaired_at = polyval(p_local_at, x(i));
                else % 窗口数据不足时使用全局拟合
                    repaired_wt = fit_wt(i);
                    repaired_at = fit_at(i);
                end
                
                % 执行修复
                cleaned_WT3(i, col) = repaired_wt;
                cleaned_AT3(i, col) = repaired_at;
            end
        end
    end
    
    % 三维可视化
    figure('Name','三维数据修复可视化');
    subplot(2,2,1);
    scatter3(WT3(:), AT3(:), 1:numel(WT3), 10, 'r','filled');
    title('原始数据分布');
    xlabel('WT3'); ylabel('AT3');
    
    subplot(2,2,2);
    scatter3(cleaned_WT3(:), cleaned_AT3(:), 1:numel(WT3), 10, 'b','filled');
    title('修复后数据分布');
    xlabel('WT3'); ylabel('AT3');
    
    subplot(2,2,[3,4]);
    plot(combined_res, 'b'); hold on;
    plot(find(outliers), combined_res(outliers), 'ro');
    title('联合残差检测结果');
    xlabel('数据索引'); ylabel('组合残差');
    legend('正常点', '异常点');
end