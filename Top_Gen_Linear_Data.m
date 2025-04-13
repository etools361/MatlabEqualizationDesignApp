% 加载数据
WT_hist = load('WT_hist');
WT_hist = WT_hist.WT_hist;
AT0_hist = load('AT0_hist');
AT0_hist = AT0_hist.AT0_hist;

% 初始化输出变量
wTabNxC1 = [];
wTabNxC3 = [];
wTabNxC5 = [];
wTabNxC7 = [];
wTabNxC9 = [];
wTabNxC11 = [];
Freq0 = logspace(log10(0.0001),log10(0.01), 10);
C0 = 1:3:7;
N0 = 1:5;
WT_new = [];
fid = fopen('eq_data_linear.m', 'w');
strName = [];
strC0 = [];
for ii=1:length(C0)
    if ii==length(C0)
        strC0  = [strC0, sprintf('%0.1f', C0(ii))];
        strName = [strName, sprintf('wTabNxC%d', C0(ii))];
    else
        strC0  = [strC0, sprintf('%0.1f,', C0(ii))];
        strName = [strName, sprintf('wTabNxC%d,', C0(ii))];
    end
end
strH = sprintf('function [CellDataAll,C0] = eq_data_linear()\n');
fprintf(fid, strH);
fprintf(fid, 'C0=[%s];\n', strC0);
for mm=1:length(C0)
    ic = 1;
    for ii=1:length(N0)
        for jj=1:ii
            [a,b,c]=size(WT_hist{mm,kk,ii});
            for kk=1:length(Freq0)
                if c>1
                    WT_new(ic+jj-1,kk) = WT_hist{mm,kk,ii}(1,1,jj);
                    WT_new(ic+ii+jj-1,kk) = AT0_hist{mm,kk,ii}(1,1,jj);
                else
                    WT_new(ic+jj-1,kk) = WT_hist{mm,kk,ii}(1,jj);
                    WT_new(ic+ii+jj-1,kk) = AT0_hist{mm,kk,ii}(1,jj);
                end
            end
        end
        ic = ic + 2*ii;
    end
    [a,b]= size(WT_new);
    fprintf(fid, '\nwTabNxC%d = {\n', C0(mm));
    ic = 1;
    for ii=1:length(N0)
        fprintf(fid, '[\n');
        for kk=1:ii
            for jj=1:b
                if jj==b
                    fprintf(fid, '%0.4f;\n', WT_new(ic,jj));
                else
                    fprintf(fid, '%0.4f,', WT_new(ic,jj));
                end
            end
            ic = ic + 1;
        end
        fprintf(fid, '];\n');
        fprintf(fid, '[\n');
        for kk=1:ii
            for jj=1:b
                if jj==b
                    fprintf(fid, '%0.4f;\n', WT_new(ic,jj));
                else
                    fprintf(fid, '%0.4f,', WT_new(ic,jj));
                end
            end
            ic = ic + 1;
        end
        fprintf(fid, '];\n');
    end
    fprintf(fid, '};');
end
fprintf(fid, '\nCellDataAll = {%s};\n', strName);
fprintf(fid, 'end\n');
fclose(fid);

