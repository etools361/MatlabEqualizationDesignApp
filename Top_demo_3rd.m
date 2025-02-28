%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 3级设计实例
%--------------------------------------------------------------------------
N = 3;
m = 0.1;
n = 10;
c = 1.5;
maxIter = 30;
tol = 1e-10;
[aSol,wSol,deltaSol,xSol] = funRemezEquirippleRational(N, m, n, c, maxIter, tol, 7);


W = wSol;
A = aSol;
syms w
n = length(A);
y = 0;
for ii=1:n
    k = sqrt((1+A(ii)/20)/(1-A(ii)/20));
    if ii==1
        a = W(ii)*k;
        b = W(ii)/k;
        y = y + 20*w^2*(a^2-b^2)/((w^2+a^2)*(w^2+b^2));
    else
        a = W(ii)*k;
        b = W(ii)/k;
        ar = W(1)^2/W(ii)*k;
        br = W(1)^2/W(ii)/k;
        y = y + 20*w^2*(a^2-b^2)/((w^2+a^2)*(w^2+b^2)) + 20*w^2*(ar^2-br^2)/((w^2+ar^2)*(w^2+br^2));
    end
end
y2 = simplify(y)
vpa(y2,3)
[n,d] = numden(y2);
zn = vpa(solve(n==0),8);
zd = vpa(solve(d==0),8);


  A_1 = 1.38852;
  A_2 = 1.51812;
  A_3 = 1.51812;
  k1 = sqrt((1+A_1/20)/(1-A_1/20));
  k2 = sqrt((1+A_2/20)/(1-A_2/20));
  k3 = sqrt((1+A_3/20)/(1-A_3/20));
  
  W_1 = 1;
  W_2 = 6.6394;
  W_3 = 1/6.6394;
  
a1 = W_1*k1;
b1 = W_1/k1;
a2 = W_2*k2;
b2 = W_2/k2;
a3 = W_3*k3;
b3 = W_3/k3;

% zobel network
ai = [a2,a1,a3];
bi = [b2,b1,b3];
wi = ai-bi;
Ki = (ai./bi);
Ci=1./(50.*wi)
Li = 50./wi
R2in1=50.*(Ki-1)./(Ki+1)
R2i=100.*Ki./(Ki.^2-1)

syms s
n1 = double(vpa(coeffs(expand((s+a2)*(s+a1)*(s+a3))),10));
d1 = double(vpa(coeffs(expand((s+b2)*(s+b1)*(s+b3))),10));
% [r,p,k]=residue(100*fliplr(n1),fliplr(d1))
[r,p,k]=residue(2*fliplr(n1)-fliplr(d1),50*fliplr(d1))
% 1/(R+s*L),1/(a*s+b)=L//R-->1/(1/R+1/(sL))=sL/(sL/R+1),L+R-->sL+R, C//R-->1/(1/R+sC) 
% 100+98.01/(x+6.15)=100+1/(x/98.1+6.15/98.1)=100+1/(0.01s+1/15.95)

1./r
r./p

r2 = r*100;

w = logspace(-2,2,1000);
% syms w real
s = 1i.*w;
H = 100./(100+1./(1/15.9+0.0102.*s)+1./(1/17.3+0.0620.*s)+1./(1/22.6+0.318.*s));
% xx=vpa(simplify(real(H*conj(H))),3)
% [n2,d2] = numden(xx)
% n3 = coeffs(n2)
% n4 = n3/n3(end)
% d3 = coeffs(d2)
% d4 = d3/d3(end)
% syms w
% vpa(solve(w^6 + n4(3)*w^4 + n4(2)*w^2 + n4(1)==0),3)
% syms w
% vpa(solve(w^6 + d4(3)*w^4 + d4(2)*w^2 + d4(1)==0),3)
semilogx(w, 20.*log10(abs(H)), '-r', 'linewidth', 2);
hold on;
semilogx([0.1,0.1], [-4,-3.4], '--g', 'linewidth', 1);
semilogx([10,10], [-4,-0.4], '--g', 'linewidth', 1);
hold off;
grid on;
xlabel('w/rad/s');
ylabel('Mag/dB');
title('Order=3, \delta=3dB, w_l=0.1 rad/s,w_h=10 rad/s');


semilogx(w(1:end-1), diff(20.*log10(abs(H)))./diff(log10(w)), '-r', 'linewidth', 2);
hold on;
semilogx(w, 1.5.*ones(size(w)), '--g', 'linewidth', 1);
hold off;
grid on;
xlabel('w/rad/s');
ylabel('dH/dB/dec');
title('Order=3, \delta=3dB, w_l=0.1 rad/s,w_h=10 rad/s');

syms Cb K w real
s = 1i*w;
H = (50*Cb*K*s-50*Cb*s+1)/(K-50*Cb*s+50*Cb*K*s);
H2 = simplify(real(conj(H)*H))
(w^2 + 1/(50^2*Cb^2*(K-1)^2))/(w^2 + K^2/(2500*Cb^2*(K-1)^2))


