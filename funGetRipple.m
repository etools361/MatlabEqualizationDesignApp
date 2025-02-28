%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 获取纹波值
%--------------------------------------------------------------------------
function rip = funGetRipple(A, m)

% A = 15;
k = sqrt((A/20+1)/(1-A/20));
% mT = linspace(1.2,10,100);
% yT = [];
% for jj=1:length(mT)
%     m = mT(jj);
    y = 0;
%     y1 = 0;
%     y2 = 0;
    for ii=-20:20
        a = m^(ii)*k;
        b = m^(ii)/k;
        yt = 20.*(a^2-b^2)./((1+a^2).*(1+b^2));
%         yt1 = 20.*(a^2)./((x.^2+a^2));
%         yt2 = 20.*(b^2)./((x.^2+b^2));
%         semilogx(x, yt, '-');
%         hold on;
        y = y + yt;
%         y1 = y1 + yt1;
%         y2 = y2 + yt2;
    end
    rip = y-20*2*log(k)/log(m);
%     yT(jj)=y;
% end
%     Am = 20*2*log(k)./log(mT);
% plot(mT, yT-Am, '--m', 'linewidth', 1);
% grid on;
% hold on;
% plot(mT, Am, '-b', 'linewidth', 2);
% semilogx(x, y1, '--g', 'linewidth', 2);
% hold off;
