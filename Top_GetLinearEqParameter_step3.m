% 差值计算
tic;

fl = 1e6;
fh = 1e9;
delta = 3;
N = 2;
[w_guess, a_guess] = calculateEQParamsLin(fl, fh, N, delta);

A0 = a_guess;
W0 = w_guess;
% x0 = logspace(log10(11.11e6),log10(1e9),10000);
x0 = linspace(fl,fh,10000);
k = fh/0.9;
[val, x] =  funEvalLinearEq(A0./k, W0.*k, x0, 0);
% hold on;
% plot(x(1:end-1), diff(val)./diff(x), '-b', 'linewidth', 2);
figure(1);
plot(x, val, '-b', 'linewidth', 2);
xlabel('Freq/Hz');
ylabel('H_{dB}/dB');
% semilogx(x(1:end-1), diff(val)./diff(x), '-b', 'linewidth', 2);
grid on;
figure(2);
semilogx(x(1:end-1), diff(val)./diff(x), '-b', 'linewidth', 2);
xlabel('Freq/Hz');
ylabel('dH/df');
grid on;

% hold off
toc;

