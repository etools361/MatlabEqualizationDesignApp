%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-01(yyyy-mm-dd)
% 使用查表插值方式计算均衡器参数
%--------------------------------------------------------------------------
[wTabNxC2, wTabNxC5, wTabNxC8, wTabNxC10] = eq_data();
R = logspace(log10(1.5),log10(10),9);
C = [2,5,8,10];
wl = 0.1/(2*pi);
wh = 10/(2*pi);
Nset = 3;
Cset = 1.5;

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
    wxc2(1,ii) = funSplineInterp(R, wTabN10C2(ii,:)', Rset);
    axc2(1,ii) = funSplineInterp(R, aTabN10C2(ii,:)', Rset);
    wxc5(1,ii) = funSplineInterp(R, wTabN10C5(ii,:)', Rset);
    axc5(1,ii) = funSplineInterp(R, aTabN10C5(ii,:)', Rset);
    wxc8(1,ii) = funSplineInterp(R, wTabN10C8(ii,:)', Rset);
    axc8(1,ii) = funSplineInterp(R, aTabN10C8(ii,:)', Rset);
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

% calculate Parameter
% Zobel-T network
k = sqrt((A./20+1)./(1-A./20));
w = W.*(k-1./k);
K = k.^2;
R1 = 50.*(K-1)./(K+1);
R2 = 100.*K./(K.^2-1);
C1 = 1./(50.*w);
L1 = 50./w;
% display
fprintf('----------Zobel Eq----------\n');
for ii=1:length(C1)
    fprintf('Stage %d:\n', ii);
    fprintf('    R%d = %s Ohm;\n', 2*ii-1, Data2Suffix(R1(ii),'0.3'));
    fprintf('    R%d = %s Ohm;\n', 2*ii, Data2Suffix(R2(ii),'0.3'));
    fprintf('    C%d = %s F;\n', ii, Data2Suffix(C1(ii),'0.3'));
    fprintf('    L%d = %s H;\n', ii, Data2Suffix(L1(ii),'0.3'));
end

% Ser RC
a = k.*W;
b = W./k;
zeros_v = -a;
poles_v = -b;
[residues_r, poles_r, constant_r] = funResidueFromRoots(zeros_v, poles_v);
residues_r2 = residues_r.*100;
Ri = residues_r2./abs(poles_r);
Ci = 1./residues_r2;
fprintf('----------Ser RC----------\n');
for ii=1:length(Ri)
    fprintf('Stage %d:\n', ii);
    fprintf('    R%d = %s Ohm;\n', ii, Data2Suffix(Ri(ii),'0.3'));
    fprintf('    C%d = %s F;\n', ii, Data2Suffix(Ci(ii),'0.3'));
end

% Pal RL
% (2*p-q)/(50*q)
residues = funComputeResiduesForRL(-a, -b);
Li  = 1./residues;
Ri2 = b'./residues;

fprintf('----------Pal RL----------\n');
for ii=1:length(Li)
    fprintf('Stage %d:\n', ii);
    fprintf('    R%d = %s Ohm;\n', ii, Data2Suffix(Ri2(ii),'0.3'));
    fprintf('    L%d = %s H;\n', ii, Data2Suffix(Li(ii),'0.3'));
end
addpath('..\MatlabFilterDesignApp')
% [strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, Apr, Asr, bw, fShape);
% [strNetlist] = funSynthesisTransAndGenNetlist2(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, cellValueNetlist);
RS = 50;
RL = 50;
strNetlistHeader = {
    'V0 V 1 0 1';
    sprintf('RS R 1 2 %f',  RS);
};
strNetlistTail = {
    sprintf('RL R %d 0 %f',Nset+2, RL);
};
strNetlistBody = {};
for ii=1:Nset
    strNetlistBody{2*ii-1,1} = sprintf('R%d R %d %d %f', ii, ii+1, ii+2, Ri(ii));
    strNetlistBody{2*ii,1}   = sprintf('C%d C %d %d %f', ii, ii+1, ii+2, Ci(ii));
end
strNetlist = [strNetlistHeader;strNetlistBody;strNetlistTail];
[iType, Value, cellNode1, CellNode2, cellName] = funSimNetlist2Array(strNetlist);
% netlist standard
[node1, node2] = funSimNetlistRenode(cellNode1, CellNode2);
% netlist analysis
[maxNode, nL, nI, nV, nR0] = funSimNetlistAna(iType, Value, node1, node2);
% 标记给定器件的非GND节点用于获取结果位置
strDevice = 'RL'; % 显示结果的器件
[retNode] = funGetDeviceNode(cellName, node1, node2, strDevice);
% generate matrix, (MM*s+MN)*MX = MV
[MM, MN, MV, MX] = funSimNetlist2Matrix(iType, Value, node1, node2, maxNode, nL, nI, nV, nR0, cellName);
f0 = 1e-2;
f1 = 1e2;
% f0 = 1e4;
% f1 = 1e6;
N     = 500;
if exist('h2','var') && ishandle(h2)
else
    h2    = figure(2);
    axMag = axes(h2);
end
if exist('h3','var') && ishandle(h3)
else
    h3    = figure(3);
    axPhase = axes(h3);
end
IdealFreq = 0;
IdealMag = 0;
IdealPhase = 0;
% [IdealFreq, IdealMag, IdealPhase, P, Z] = funSimFilterIdeal(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape, f0, f1, N);
funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, IdealFreq, IdealMag, IdealPhase, 1);
ylim(axMag, [-12,-6])
% axMag = [];
% axPhase = [];
% funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, IdealFreq, IdealMag, IdealPhase, logscaleEn);
