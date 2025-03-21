N = 4;
m = 0.1;
n = 10;
delta = 3;
c = delta/log10(n/m);
maxIter = 20;
tol = 1e-6;
[aSol,wSol,deltaSol,xSol] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, 0);

W = wSol;
A = aSol;
syms w
n = floor((N+1)/2);
y = 0;
for ii=1:n
    k = sqrt((1+A(ii)/20)/(1-A(ii)/20));
    if mod(N,2) && ii==1
        a = W(ii)*k;
        b = W(ii)/k;
        y = y + 20*w^2*(a^2-b^2)/((w^2+a^2)*(w^2+b^2));
    else
        a = W(ii)*k;
        b = W(ii)/k;
%         if ~mod(N,2) && ii==1
            ar = WC^2/W(ii)*k;
            br = WC^2/W(ii)/k;
%         else
%             ar = W(1)^2/W(ii)*k;
%             br = W(1)^2/W(ii)/k;
%         end
        y = y + 20*w^2*(a^2-b^2)/((w^2+a^2)*(w^2+b^2)) + 20*w^2*(ar^2-br^2)/((w^2+ar^2)*(w^2+br^2));
    end
end
%
w0 = logspace(-1,1,100);
yreal = [];
for ii=1:length(w0)
    yreal(ii) = double(subs(y,w,w0(ii)));
end
semilogx(w0, yreal, '-r', 'linewidth', 2);
grid on;
hold off;

y2    = simplify(y);
y2_jw0 = subs(y2,w,1i*w);
y2_jw  = simplify(y2_jw0);
dy2   = simplify(diff(y2));
dy2_jw0 = subs(dy2,w,1i*w);
dy2_jw   = simplify(dy2_jw0);

[n,d] = numden(dy2);
zn = vpa(solve(n==0),4)
zd = vpa(solve(d==0),4)

plot((real(zn)),(imag(zn)),'ob');
grid on;
hold on;
plot((real(zd)),(imag(zd)),'xr');
hold off;
axis equal
xlabel('imag');
ylabel('real');


%----------------------------------
W = logspace(-1,1,120);
for jj=1:length(zn)
    if jj==1
        ZN3 = double(W+zn(jj));
    else
        ZN3 = double(ZN3.*(W+zn(jj)));
    end
end
for jj=1:length(zd)
    if jj==1
        ZD3 = double((W+zd(jj)));
    else
        ZD3 = double(ZD3.*(W+zd(jj)));
    end
end
YY3 = ZN3./ZD3*1e3;
plot3(W,zeros(size(W)),log10(abs(YY3)),'-r');
grid on;
%-------------------------
x = linspace(-1e1,1e1,120);
y = linspace(0,1e1,60);
[X,Y] = meshgrid(x,y);
W = X+1i.*Y;
for jj=1:length(zn)
    if jj==1
        ZN = double(W+zn(jj));
    else
        ZN = double(ZN.*(W+zn(jj)));
    end
end
for jj=1:length(zd)
    if jj==1
        ZD = double((W+zd(jj)));
    else
        ZD = double(ZD.*(W+zd(jj)));
    end
end
YY = ZN./ZD*1e3;
aY = abs(YY);
surf(X,Y,log10(aY));
set(gcf,'color', [1,1,1]);
grid off;
% plot((real(zn)),(imag(zn)),'ob');
% grid on;
% hold on;
% plot((real(zd)),(imag(zd)),'xr');
% hold off;
axis equal



P1 = [];
NP1 = 1:mNN;
for jj=1:length(R)
    for ii=1:1:mNN
        [w0, iw0] = sort(abs(wSol{jj,ii}));
%         AA = abs(aSol{jj,ii});
        AA = w0;
        AA0 = AA(end);
        if ii==1
            A00 = AA0;
        end
        P1(jj,ii) = log10(AA0./A00);
    end
end
NN = NN(2:end);
P1(:,1)=[];
[N3, R3] = meshgrid(NN,R);
surf(N3, R3, P1);



