function [val, x] = funEvalLinearEq(A, W, x, DispEn)
if nargin<3
    x = linspace(0,1,1000);
else
    if isempty(x)
        x = linspace(0,1,1000);
    end
end
if nargin<4
    DispEn = 1;
end
xx = x.^2;
val = 0;
aa = log(10);
a = W .* sqrt((30./aa+A.*W)./(10./aa-A.*W));
b = W .* sqrt((30./aa-A.*W)./(10./aa+A.*W));
for i = 1:length(a)
    val = val + 10 .* log10((xx+b(i)^2) ./ (xx+a(i)^2));
end
if DispEn
    plot(x, val, '-r', 'linewidth', 2);
    grid on;
end
end