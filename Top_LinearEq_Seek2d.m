N = 3;      % 项数
% m = 0.025;      % 区间左端点
% n = 0.55;      % 区间右端点
cT = 3:1:10;    % 目标中心值
mCT = length(cT);
AT1 = [];
AT2 = [];
AT3 = [];
WT1 = [];
WT2 = [];
WT3 = [];
ATHist = [];
WTHist = [];
deltaTHist = [];
for kk=1:mCT
    c = cT(kk);
    tol = 1e-6; % 收敛容差
    max_iter = 50; % 最大迭代次数

    % 调用函数
    % [delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter);
    kT = linspace(1.5,10,100);
    nn = length(kT);
    deltaT=[];
    AT    =[];
    WT    =[];
    for ii=1:nn
        if isempty(ATHist)
            A_set = A;
            W_set = W;
            delta_set = delta;
        else
            A_set = ATHist(ii,:);
            W_set = WTHist(ii,:);
            delta_set = deltaTHist(ii);
        end
%         if ~mod(kk,2)
            k = kT(ii);
%         else
%             k = kT(nn-ii+1);
%         end
        p = 1/k*0.9;
        m = 1/k*p;    % 区间左端点
        n = k*p;      % 区间右端点
        [delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, A_set, W_set,delta_set, 0);
        deltaT(ii) = delta;
        AT(ii,:)   = A;
        WT(ii,:)   = W;
        plot(kT(1:ii), deltaT(1:ii), '-r', 'linewidth', 2);
        grid on;
        drawnow;
    end
    ATHist = AT;
    WTHist = WT;
    deltaTHist = deltaT;
    AT1(kk,:) = AT(:,1);
    AT2(kk,:) = AT(:,2);
    AT3(kk,:) = AT(:,3);
    WT1(kk,:) = WT(:,1);
    WT2(kk,:) = WT(:,2);
    WT3(kk,:) = WT(:,3);
end
% AT3
% p1 =        6.79;
% p2 =       48.02;
% p3 =      -35.46;
% q1 =       12.49;
% q2 =      -9.164;

