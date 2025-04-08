% 修复系数


WT_hist  = load('WT0');
WT_hist = WT_hist.WT0;
AT0_hist = load('AT0');
AT0_hist = AT0_hist.AT0;
% deltaT0 = deltaT0;%(C0,Freq,N)

Freq0 = logspace(log10(0.001),log10(0.2), 10);
Freq1 = 0.9;
C0 = 1:2:11;
N0 = 1:5;
for jj=1:length(N0)
    N = N0(jj);
    WT_new = [];
    AT_new = [];
    for ll=1:N
        for kk=1:length(N0)
            for ii=1:length(C0)
                WT_new(ii,kk,ll) = WT_hist{ii,kk,N}(ll);
                AT_new(ii,kk,ll) = AT0_hist{ii,kk,N}(ll);
            end
        end
    end

    % fix bugs
    [a,b,c] = size(WT_new);
    WT2_repaired = [];
    AT2_repaired = [];
    for ii=1:N
        WT2 = WT_new(:,:,ii);
        AT2 = AT_new(:,:,ii);
        [WT2_repaired(:,:,ii), AT2_repaired(:,:,ii)] = auto_repair_data(WT2, AT2);
    end

    for ll=1:length(C0)
        for kk=1:length(N0)
            WT_hist{ll,kk,N} = WT2_repaired(ll,kk,:);
            AT0_hist{ll,kk,N} = AT2_repaired(ll,kk,:);
        end
    end

end

iii = 4;
% semilogx(WT_new(:,:,iii), AT_new(:,:,iii),'*')
grid on;
% hold on;
semilogx(WT2_repaired(:,:,iii), AT2_repaired(:,:,iii),'o')
% hold off;

% tic;
% m = 0.001*2;% 区间左端点
% n = 0.9*2;
% c = 2;    % 目标中心值
% N = 5;
% tol = 1e-10; % 收敛容差
% max_iter = 100; % 最大迭代次数
% % 计算范围
% ic  = interp1(C0, 1:length(C0),c);
% iff = funinterp1(Freq0, 1:length(Freq0),m/n*0.9);
% a_guess = [];w_guess = [];
% for ii=1:N
%     a_guess(ii) = AT0_hist{ic,iff,N}(ii);
%     w_guess(ii) = WT_hist{ic,iff,N}(ii);
% end
% delta_guess = 0.1;
% dispEn = 3;
% [delta0, A0, W0, x_points0] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, a_guess, w_guess, delta_guess, dispEn);
% toc;
save('AT0_hist','AT0_hist');
save('WT_hist','WT_hist');


