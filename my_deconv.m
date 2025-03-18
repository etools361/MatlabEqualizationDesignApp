function [quotient, remainder] = my_deconv(num, den)
    % 输入：num 是分子多项式系数，den 是分母多项式系数
    % 输出：quotient 是商，remainder 是余数

    % 初始化
    quotient = [];
    remainder = num;

    % 长除法
    while length(remainder) >= length(den)
        % 计算当前项的系数
        current_term = remainder(1) / den(1);
        quotient = [quotient, current_term]; % 添加到商中

        % 更新余数
        remainder = remainder - my_conv(den, [current_term, zeros(1, length(remainder) - length(den))]);
        remainder = remainder(2:end); % 去掉最高次项
    end
end