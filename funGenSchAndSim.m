%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2025-03-18(yyyy-mm-dd)
% 由A和W值综合出原理图值并执行仿真
% Slope = 0;% 0:Positive Slope;1:Negative Slope
% Type = 1;% netlist:0,zobel network, 1,RC/RL Serial, 2, RC/RL Parallel
% f0 = 1e-2;
% f1 = 1e2;
%--------------------------------------------------------------------------
function funGenSchAndSim(A, W, Rs, Slope, Type, fl, fu)
global h2;
global h3;
% calculate Parameter
% Zobel-T network
Nset = length(A);
k = sqrt((A./20+1)./(1-A./20));
% a = k.*W;% k = sqrt(a/b), k' = sqrt(b/a)
% b = W./k;% W = sqrt(a*b), W' = sqrt(a*b)
if Slope
    a  = k.*W;
    K  = k.^2;
    w0 = a/(K-1);
    R1 = Rs.*(K-1)./(K+1);
    R2 = 2*Rs.*K./(K.^2-1);
    C1 = 1./(Rs.*w0);
    L1 = Rs./w0;
else
    w  = W.*(k-1./k);
    K  = k.^2;
    R1 = Rs.*(K-1)./(K+1);
    R2 = 2*Rs.*K./(K.^2-1);
    C1 = 1./(Rs.*w);
    L1 = Rs./w;
end
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
if Slope
    zeros_v  = -b;
    poles_v  = -a;
    Qs = my_poly(-a);
    Ps = my_poly(-b);
    Q0 = my_prod(a);
    P0 = my_prod(b);
    PP = Ps.*Q0./P0-Qs;
%     [residues_r, poles_r] = my_residue(zeros_v2(1:end-1), -a);
    [residues_r, poles_r]=my_residue((PP(1:end-1)),-a);
    residues_r2 = residues_r.*2.*Rs;
    Ri = residues_r2;
    Ci = residues_r2./a;% Li
else
    zeros_v = -a;
    poles_v = -b;
    zeros_v2 = my_poly(zeros_v);
    [residues_r, poles_r]=my_residue(zeros_v2,-b);
%     [residues_r, poles_r, constant_r] = funResidueFromRoots(zeros_v, poles_v);
    residues_r2 = residues_r.*2.*Rs;
    Ri = residues_r2./abs(poles_r);
    Ci = 1./residues_r2;
end
fprintf('----------Ser RC----------\n');
for ii=1:length(Ri)
    fprintf('Stage %d:\n', ii);
    fprintf('    R%d = %s Ohm;\n', ii, Data2Suffix(Ri(ii),'0.3'));
    fprintf('    C%d = %s F;\n', ii, Data2Suffix(Ci(ii),'0.3'));
end

% Pal RL
% (2*p-q)/(50*q)
if Slope
    Ps = my_poly((-b));
    Qs = my_poly((-a));
    P0 = my_prod(b);
    Q0 = my_prod(a);
    QQ = Qs;
    PP = Ps*Q0./P0-Qs;
    [residues_r, poles_r] = my_residue(PP(1:end-1), -a);
%     [residues_r,poles_r,constant_r]=residue(PP(1:end-1)./PP0,QQ);
    residues_r2 = residues_r.*2./Rs;
    Ri2 = 1./(residues_r2);
    Li  = -1./(poles_r.*Ri2);% Ci
else
%     residues2 = funComputeResiduesForRL(-a, -b, Rs);
    Ps = my_poly((-b));
    Qs = my_poly((-a));
    PP = Qs-Ps;
    [residues, poles_r] = my_residue(PP(2:end), -b);
    residues2 = residues.*2./Rs;
    Li  = 1./residues2;
    Ri2 = b./residues2;
end

fprintf('----------Pal RL----------\n');
for ii=1:length(Li)
    fprintf('Stage %d:\n', ii);
    fprintf('    R%d = %s Ohm;\n', ii, Data2Suffix(Ri2(ii),'0.3'));
    fprintf('    L%d = %s H;\n', ii, Data2Suffix(Li(ii),'0.3'));
end
addpath('..\MatlabFilterDesignApp')
% [strNetlist] = funSynthesisFilter(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, Apr, Asr, bw, fShape);
% [strNetlist] = funSynthesisTransAndGenNetlist2(fType, fShape, TeeEn, n, Rs, Rl, fp, bw, cellValueNetlist);
RS = Rs;
RL = Rs;
strNetlistHeader = {
    'V0 V 1 0 1';
    sprintf('RS R 1 2 %f',  RS);
};
strNetlistBody = {};
if Type==0% zobel network
    strNetlistTail = {
        sprintf('RL R %d 0 %f',3*Nset+2, RL);
    };
    if Slope
        for ii=1:Nset
            strNetlistBody{5*(ii-1)+1,1} = sprintf('R%d R %d %d %f', ii,    3*(ii-1)+2, 3*(ii-1)+3, R1(ii));
            strNetlistBody{5*(ii-1)+2,1}   = sprintf('R%d R %d %d %f', ii+10, 3*(ii-1)+3, 3*(ii-1)+5, R1(ii));
            strNetlistBody{5*(ii-1)+3,1}   = sprintf('R%d R %d %d %f', ii+20, 3*(ii-1)+3, 3*(ii-1)+4, R2(ii));
            strNetlistBody{5*(ii-1)+4,1}   = sprintf('L%d L %d %d %f', ii,    3*(ii-1)+2, 3*(ii-1)+5, L1(ii));
            strNetlistBody{5*(ii-1)+5,1}   = sprintf('C%d C %d %d %f', ii,    3*(ii-1)+4, 0, C1(ii));
        end
    else
        for ii=1:Nset
            strNetlistBody{5*(ii-1)+1,1} = sprintf('R%d R %d %d %f', ii,    3*(ii-1)+2, 3*(ii-1)+3, R1(ii));
            strNetlistBody{5*(ii-1)+2,1}   = sprintf('R%d R %d %d %f', ii+10, 3*(ii-1)+3, 3*(ii-1)+5, R1(ii));
            strNetlistBody{5*(ii-1)+3,1}   = sprintf('R%d R %d %d %f', ii+20, 3*(ii-1)+3, 3*(ii-1)+4, R2(ii));
            strNetlistBody{5*(ii-1)+4,1}   = sprintf('C%d C %d %d %f', ii,    3*(ii-1)+2, 3*(ii-1)+5, C1(ii));
            strNetlistBody{5*(ii-1)+5,1}   = sprintf('L%d L %d %d %f', ii,    3*(ii-1)+4, 0, L1(ii));
        end
    end
elseif Type==1 % RC/RL network
    strNetlistTail = {
        sprintf('RL R %d 0 %f',Nset+2, RL);
    };
    if Slope
        for ii=1:Nset
            strNetlistBody{2*ii-1,1} = sprintf('R%d R %d %d %f', ii, ii+1, ii+2, Ri(ii));
            strNetlistBody{2*ii,1}   = sprintf('L%d L %d %d %f', ii, ii+1, ii+2, Ci(ii));
        end
    else
        for ii=1:Nset
            strNetlistBody{2*ii-1,1} = sprintf('R%d R %d %d %f', ii, ii+1, ii+2, Ri(ii));
            strNetlistBody{2*ii,1}   = sprintf('C%d C %d %d %f', ii, ii+1, ii+2, Ci(ii));
        end
    end
else
    strNetlistTail = {
        sprintf('RL R %d 0 %f',2, RL);
    };
    if Slope
        for ii=1:Nset
            strNetlistBody{2*ii-1,1} = sprintf('R%d R %d %d %f', ii, 2, ii+2, Ri2(ii));
            strNetlistBody{2*ii,1}   = sprintf('C%d C %d %d %f', ii, ii+2, 0, Li(ii));
        end
    else
        for ii=1:Nset
            strNetlistBody{2*ii-1,1} = sprintf('R%d R %d %d %f', ii, 2, ii+2, Ri2(ii));
            strNetlistBody{2*ii,1}   = sprintf('L%d L %d %d %f', ii, ii+2, 0, Li(ii));
        end
    end
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
% f0 = 1e4;
% f1 = 1e6;
N     = 500;
if exist('h2','var') && ~~(ishandle(h2))
    axMag = axes(h2);
else
    h2    = figure(2);
    axMag = axes(h2);
end
if exist('h3','var') && ~~(ishandle(h3))
    axPhase = axes(h3);
else
    h3    = figure(3);
    axPhase = axes(h3);
end
IdealFreq = 0;
IdealMag = 0;
IdealPhase = 0;
% [IdealFreq, IdealMag, IdealPhase, P, Z] = funSimFilterIdeal(fType, TeeEn, n, Rs, Rl, fp, fs, Ap, As, bw, fShape, f0, f1, N);
funACSim(axMag, axPhase, fl, fu, N, cellName, MX, MM, MN, MV, Value, node1, node2, IdealFreq, IdealMag, IdealPhase, 1);
ylim(axMag, [-12,-6])
% axMag = [];
% axPhase = [];
% funACSim(axMag, axPhase, f0, f1, N, cellName, MX, MM, MN, MV, Value, node1, node2, IdealFreq, IdealMag, IdealPhase, logscaleEn);

end