function [w_guess, a_guess] = calculateEQParamsLin(fl, fh, N, delta)
WT_hist  = load('WT_hist');
WT_hist  = WT_hist.WT_hist;
AT0_hist = load('AT0_hist');
AT0_hist = AT0_hist.AT0_hist;

f1 = 0.9*fl/fh;
m = f1; % 区间左端点
n = 0.9;
c = delta/(n-m);       % 目标中心值
f_target = m/n * 0.9;
% 定义已知的校准点
Freq0 = logspace(log10(0.0001),log10(0.01), 10);
C0 = [1,4,7];
C0 = C0(:)';
Freq0 = Freq0(:)';
if c <= C0(1)
    ic_low = 1;
    ic_high = 2;
    alpha_c = (c - C0(1))/(C0(2) - C0(1));
elseif c >= C0(end)
    ic_low = length(C0)-1;
    ic_high = length(C0);
    alpha_c = (c - C0(end-1))/(C0(end) - C0(end-1));
else
    ic_low = find(C0 <= c, 1, 'last');
    ic_high = ic_low + 1;
    alpha_c = (c - C0(ic_low))/(C0(ic_high) - C0(ic_low));
end

% 处理频率的插值
% 找到Freq0中最近的上下界索引
if f_target <= Freq0(1)
    iff_low = 1;
    iff_high = 2;
    alpha_f = (f_target - Freq0(1))/(Freq0(2) - Freq0(1));
elseif f_target >= Freq0(end)
    iff_low = length(Freq0)-1;
    iff_high = length(Freq0);
    alpha_f = (f_target - Freq0(end-1))/(Freq0(end) - Freq0(end-1));
else
    iff_low = find(Freq0 <= f_target, 1, 'last');
    iff_high = iff_low + 1;
    alpha_f = (f_target - Freq0(iff_low))/(Freq0(iff_high) - Freq0(iff_low));
end

% 双线性插值获取初始猜测值
a_guess = zeros(1,N);
w_guess = zeros(1,N);
for ii = 1:N
    % 获取四个角点的值
    a11 = AT0_hist{ic_low, iff_low, N}(ii);
    a12 = AT0_hist{ic_low, iff_high, N}(ii);
    a21 = AT0_hist{ic_high, iff_low, N}(ii);
    a22 = AT0_hist{ic_high, iff_high, N}(ii);
    
    w11 = WT_hist{ic_low, iff_low, N}(ii);
    w12 = WT_hist{ic_low, iff_high, N}(ii);
    w21 = WT_hist{ic_high, iff_low, N}(ii);
    w22 = WT_hist{ic_high, iff_high, N}(ii);
    
    % 进行双线性插值
    a_guess(ii) = (1-alpha_c)*(1-alpha_f)*a11 + ...
                 (1-alpha_c)*alpha_f*a12 + ...
                 alpha_c*(1-alpha_f)*a21 + ...
                 alpha_c*alpha_f*a22;
    
    w_guess(ii) = (1-alpha_c)*(1-alpha_f)*w11 + ...
                 (1-alpha_c)*alpha_f*w12 + ...
                 alpha_c*(1-alpha_f)*w21 + ...
                 alpha_c*alpha_f*w22;
end
end