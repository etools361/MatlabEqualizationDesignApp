function [x,y] = funLogInterp1(x0, y0)
a = 1.1;
x1 = log(x0(1));
x2 = log(x0(2));
y1 = y0(1);
y2 = y0(2);
x3 = (x2-x1+a.*x2)./a;
y3 = (y2-y1+a.*y2)./a;
x  = exp(x3);
y  = y3;
end