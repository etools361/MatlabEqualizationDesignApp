%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 使用查表插值方式计算均衡器参数
%--------------------------------------------------------------------------
[wTabNxC2, wTabNxC5, wTabNxC8, wTabNxC10] = eq_data();
R = logspace(log10(1.5),log10(10),9);
C = [2,5,8,10];
Rs = 50;
Rl = Rs;
wl = 0.1*(2*pi);
wh = 10*(2*pi);
Nset = 3;
Cset = 1.5;
Slope = 0;% 0:Positive Slope;1:Negative Slope
Type = 1;% netlist:0,zobel network, 1,RC/RL Serial, 2, RC/RL Parallel

Rset = log10(wh/wl);
wTabN10C2 = wTabNxC2{Nset*2-1};
aTabN10C2 = wTabNxC2{Nset*2};
wTabN10C5 = wTabNxC5{Nset*2-1};
aTabN10C5 = wTabNxC5{Nset*2};
wTabN10C8 = wTabNxC8{Nset*2-1};
aTabN10C8 = wTabNxC8{Nset*2};
wTabN10C10 = wTabNxC10{Nset*2-1};
aTabN10C10 = wTabNxC10{Nset*2};
wx = [];wxc2=[];wxc5=[];wxc8=[];wxc10=[];
ax = [];axc2=[];axc5=[];axc8=[];axc10=[];
[a1,b1] = size(wTabN10C2);
for ii=1:a1
    wxc2(1,ii)  = funSplineInterp(R, wTabN10C2(ii,:)', Rset);
    axc2(1,ii)  = funSplineInterp(R, aTabN10C2(ii,:)', Rset);
    wxc5(1,ii)  = funSplineInterp(R, wTabN10C5(ii,:)', Rset);
    axc5(1,ii)  = funSplineInterp(R, aTabN10C5(ii,:)', Rset);
    wxc8(1,ii)  = funSplineInterp(R, wTabN10C8(ii,:)', Rset);
    axc8(1,ii)  = funSplineInterp(R, aTabN10C8(ii,:)', Rset);
    wxc10(1,ii) = funSplineInterp(R, wTabN10C10(ii,:)', Rset);
    axc10(1,ii) = funSplineInterp(R, aTabN10C10(ii,:)', Rset);
end
% [fitresult, gof] = createFitForW(R, wTabN10C2(ii,:));
wxc0 = [];
axc0 = [];
for ii=1:a1
    wxc0(ii) = funSplineInterp(C, [wxc2(ii),wxc5(ii),wxc8(ii),wxc10(ii)], Cset);
    axc0(ii) = funSplineInterp(C, [axc2(ii),axc5(ii),axc8(ii),axc10(ii)], Cset);
end

n1 = wh;
m1 = wl;
Wc = sqrt(wh*wl);
xi = logspace(log10(m1/1), log10(n1*1), 10000);
[A, W, delta] = funGetFullPara([axc0';Wc.*10.^(wxc0'*Rset/2);0], sqrt(m1*n1), Nset);
yVal = funCalcY(A, W, xi);
semilogx(xi, yVal, '-', 'linewidth', 2);
grid on;
funGenSchAndSim(A, W, Rs, Slope, Type, wl/2/pi, wh/2/pi);
