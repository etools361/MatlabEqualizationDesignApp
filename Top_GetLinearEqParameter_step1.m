% 生成系数


% m = 0.001;      % 区间左端点
% c = 1;    % 目标中心值
tol = 1e-10; % 收敛容差
max_iter = 100; % 最大迭代次数
tic;

% [delta, A, W, x_points] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter);

AT0 = [];
WT0 = [];
deltaT0 = [];
ASetT = [];
WSetT = [];
Freq = logspace(log10(0.001),log10(0.2), 10);
pp = 1000;
C0 = 1:2:11;
N0 = 1:5;
mC = length(C0);
mFreq = length(Freq);
for ll=1:mC
    c = C0(ll);
    for kk=1:mFreq
        m = Freq(kk); 
        n = 0.9;      % 区间右端点
        for ii=1:length(N0)
            N = ii;
            if ii<4
                [delta0, A0, W0, x_points0] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter);
                deltaT0(ll,kk,ii) = delta0;
                AT0{ll,kk,ii}     = A0;
                WT0{ll,kk,ii}     = W0;
            else
                nIn = floor((N-1)/2);
                for jj=1:nIn
                    xPio1 = [];xPie1 = [];yPio1 = [];yPie1 = [];
                    nIn2 = N-2*jj+1;
                    for kk2=1:nIn2
                        Nx = N-kk2;
                        xP1   = WT0{ll,kk,Nx};
                        yP1   = AT0{ll,kk,Nx};
                        xPio1(kk2) = xP1(jj); 
                        xPie1(kk2) = xP1(Nx-jj+1);
                        yPio1(kk2) = yP1(jj);
                        yPie1(kk2) = yP1(Nx-jj+1);
                        if kk2>1
                            break;
                        end
                    end
                    nn = length(xPio1);
                    [WT0{ll,kk,ii}(jj),AT0{ll,kk,ii}(jj)] = funLogInterp1(fliplr(xPio1), fliplr(yPio1));
                    [WT0{ll,kk,ii}(N-jj+1),AT0{ll,kk,ii}(N-jj+1)] = funLogInterp1(fliplr(xPie1), fliplr(yPie1));
        %             tn = 1:nn;
        %             [p0] = polyfit(tn,log(xPio1),pN);
        %             [p1] = polyfit(tn,yPio1,pN);
        %             [p2] = polyfit(tn,log(xPie1),pN);
        %             [p3] = polyfit(tn,yPie1,pN);
        %             n_next = 0; 
        %             WT0{ll,kk,ii}(N-jj+1) = exp(polyval(p2,n_next));
        %             WT0{ll,kk,ii}(jj) = exp(polyval(p0,n_next));
        %             n_next = 0; 
        %             AT0{ll,kk,ii}(jj) = polyval(p1,n_next);
        %             AT0{ll,kk,ii}(N-jj+1) = polyval(p3,n_next);
        %             semilogx([xPio1,xPio2],[yPio1,yPio2], '*');grid on;
        %             hold on;
        %             semilogx([xPie1,xPie2],[yPie1,yPie2], '*');grid on;
        %             semilogx([xo],[yo], '*');grid on;
        %             semilogx([xe],[ye], '*');grid on;
        %             hold off;
                end
                if mod(N,2)
                    % pred. 1 point
                    x1 = exp((log(WT0{ll,kk,ii}(jj))+log(WT0{ll,kk,ii}(N-jj+1)))/2);
                    y1 = (AT0{ll,kk,ii}(jj)+AT0{ll,kk,ii}(N-jj+1))/2;
                    WT0{ll,kk,ii}(nIn+1) = x1;
                    AT0{ll,kk,ii}(nIn+1) = y1;
                else
                    % pred. 2 points
                    x1 = exp((2*log(WT0{ll,kk,ii}(jj))+log(WT0{ll,kk,ii}(N-jj+1)))/3);
                    y1 = (2*AT0{ll,kk,ii}(jj)+AT0{ll,kk,ii}(N-jj+1))/3;
                    WT0{ll,kk,ii}(nIn+1) = x1;
                    AT0{ll,kk,ii}(nIn+1) = y1;
                    x1 = exp((log(WT0{ll,kk,ii}(jj))+2*log(WT0{ll,kk,ii}(N-jj+1)))/3);
                    y1 = (AT0{ll,kk,ii}(jj)+2*AT0{ll,kk,ii}(N-jj+1))/3;
                    WT0{ll,kk,ii}(nIn+2) = x1;
                    AT0{ll,kk,ii}(nIn+2) = y1;
                end
                A_set = AT0{ll,kk,ii};
                W_set = WT0{ll,kk,ii};
                ASetT{ii} = A_set;
                WSetT{ii} = W_set;
                [x0,y0] = funLogInterp1([deltaT0(ll,kk,ii-2),deltaT0(ll,kk,ii-1)], [ii-2,ii-1]);
                delta_set = x0;
                for mmmm=1:10
                    [delta, A, W, x_points, Err] = funRemezEquirippleRational_linefreq(N, m, n, c, tol, max_iter, A_set, W_set,delta_set, 3);
                    if ~Err
                        break;
                    else
                        A_set = randi(c*2,1,5);
                        if min(W_set)<0 || max(W_set)>2*n
                            W_set = logspace(log10(m), log10(n), N);
                        end
                    end
                end
        %         semilogx(W,A, '-*');grid on;
                AT0{ll,kk,ii} = A;
                WT0{ll,kk,ii} = W;
                deltaT0(ll,kk,ii) = delta;
            end
        end
    end
end
semilogx(WT0{1,1,1},AT0{1,1,1}, '*');grid on;
x1 = [];
y1 = [];
x1(1) = WT0{1,1}(1);
y1(1) = AT0{1,1}(1);
hold on;
% WT4 = WT0;
% AT4 = AT0;
for ll=1:mC%mC
    for kk=1:3%mFreq
        for ii=3%N=1.2,3,4,5
            if max(AT0{ll,kk,ii})>20
                continue;
            end
            semilogx(WT0{ll,kk,ii},AT0{ll,kk,ii}, '*');grid on;
        %     x1(ii) = WT0{ll,kk,ii}(1);
        %     y1(ii) = AT0{ll,kk,ii}(1);
        %     semilogx(WSetT{ii},ASetT{ii}, '-o');grid on;
        end
    end
end
% x1_log = log(x1);
% m = 1:N;
hold off
toc;

% WT_hist = WT0;
% AT0_hist = AT0;
% deltaT0 = deltaT0;
save('AT0','AT0');
save('WT0','WT0');


