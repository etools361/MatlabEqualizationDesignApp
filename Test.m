syms c delta m n a b w A W real
% y = 20*w*(a^2-b^2)/((w^2+a^2)*(w^2+b^2));
% ym = subs(y,w,m);
% yn = subs(y,w,n);
% a0 = W*sqrt((30+A*W)/(10-A*W));
% b0 = W*sqrt((30-A*W)/(10+A*W));
% y1 = subs(y, a, a0);
% y2 = subs(y1,b, b0);
% y3 = simplify(y2)
% pretty(y3)
y = 20*w*((10 + A*W)/((10 + A*W)*w^2 + (30 - A*W)*W^2) - (10 - A*W)/((10 - A*W)*w^2 + (30 + A*W)*W^2))
ym = subs(y,w,m);
yn = subs(y,w,n);
r = solve(ym==yn,A)
% A =  (10*((- W^4 - 2*W^2*m*n + m^3*n + m^2*n^2 + m*n^3)*(- 9*W^4 + 6*W^2*m*n + m^3*n + m^2*n^2 + m*n^3))^(1/2))/(- W^5 - 2*W^3*m*n + W*m^3*n + W*m^2*n^2 + W*m*n^3)
yy = ym+A-2*c
[n,d] = numden(yy);
y3 = collect(n,A)
simplify(n)
% r2 = solve(n==0,A)
% A = 

y3 = A^3 - 2*c*A^2 - 100*( 9*W^2 +16*W*m + 6*m^2 + m^4/W^2)/(W^2 -  m^2)^2*A +  200*c*(3*W^2 + m^2)^2/(W^2*(W^2 -  m^2)^2)
A = g(W,m,c)

a = 2;
b = 1;
x = linspace(0,10,100);
y = 20.*x.*(a^2-b^2)./((x.^2+a^2).*(x.^2+b^2))
plot(x, y, '-r', 'linewidth', 2);
grid on;

A = 10.*(a^2+b^2+sqrt(a^4+14*a^2*b^2+b^4))./(a*(a^2-b^2))

syms a b
% r = solve(a^2==W^2*(30+A*W)/(10-A*W),b^2==W^2*(30-A*W)/(10+A*W), A, W)

% A2     = simplify(A)
% [n, d] = numden(A2)
% n2     = simplify(n)
KK = sqrt(a^4+14*a^2*b^2+b^4);
% y = 5*sqrt(6)*sqrt(KK-(a^2+b^2))*(a^4+b^4-10*a^2*b^2+KK*(a^2+b^2))/(6*a^2*b^2*(a^2-b^2))

% syms a b m n

% r = solve(20*m*(a^2-b^2)/((m^2+a^2)*(m^2+b^2))==20*n*(a^2-b^2)/((n^2+a^2)*(n^2+b^2)), a)
% a = (-(m*n*(b^2 + m^2 + m*n + n^2))/(- b^2 + m*n))^(1/2);
% a^2 = m*n*(b^2 + m^2 + m*n + n^2)/(b^2 - m*n)
% y1 = 20*m*(m*n*(b^2 + m^2 + m*n + n^2)/(b^2 - m*n)-b^2)/((m^2+m*n*(b^2 + m^2 + m*n + n^2)/(b^2 - m*n))*(m^2+b^2))

