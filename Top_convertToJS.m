function Top_convertToJS()
[CellDataAll, C0] = eq_data_linear();

% 打开文件准备写入
fid = fopen('eq_data_linear.js', 'w');
if fid == -1
    error('无法创建文件 eq_data.js');
end

% 写入 JavaScript 文件头部
fprintf(fid, '// eq_data_linear.js\n');
fprintf(fid, '// 自动生成的 JavaScript 数据文件\n\n');
fprintf(fid, 'module.exports = {\n');
% C0 = [2,5,8,10];
% C0 = 1:2:11;
% 将 MATLAB 数组转换为 JavaScript 格式
for ii=1:length(C0)
    fprintf(fid, '  wTabNxC%d: [\n', C0(ii));
    convertCellArrayToJS(fid, CellDataAll{ii});
    fprintf(fid, '  ],\n\n');
end

% 写入 JavaScript 文件尾部
fprintf(fid, '};\n');

% 关闭文件
fclose(fid);
disp('eq_data.js 文件已生成');

% 辅助函数：将 MATLAB 的 cell 数组转换为 JavaScript 格式
function convertCellArrayToJS(fid, cellArray)
    for i = 1:length(cellArray)
        fprintf(fid, '    [\n');
        row = cellArray{i};
        for j = 1:size(row, 1)
            fprintf(fid, '      [');
            fprintf(fid, '%.3f,', row(j, 1:end-1));
            fprintf(fid, '%.3f', row(j, end));
            fprintf(fid, '],\n');
        end
        fprintf(fid, '    ],\n');
    end
