function Top_convertToJS()
[wTabNxC2, wTabNxC5, wTabNxC8, wTabNxC10] = eq_data();

% 打开文件准备写入
fid = fopen('eq_data.js', 'w');
if fid == -1
    error('无法创建文件 eq_data.js');
end

% 写入 JavaScript 文件头部
fprintf(fid, '// eq_data.js\n');
fprintf(fid, '// 自动生成的 JavaScript 数据文件\n\n');
fprintf(fid, 'module.exports = {\n');

% 将 MATLAB 数组转换为 JavaScript 格式
fprintf(fid, '  wTabNxC2: [\n');
convertCellArrayToJS(fid, wTabNxC2);
fprintf(fid, '  ],\n\n');

fprintf(fid, '  wTabNxC5: [\n');
convertCellArrayToJS(fid, wTabNxC5);
fprintf(fid, '  ],\n\n');

fprintf(fid, '  wTabNxC8: [\n');
convertCellArrayToJS(fid, wTabNxC8);
fprintf(fid, '  ],\n\n');

fprintf(fid, '  wTabNxC10: [\n');
convertCellArrayToJS(fid, wTabNxC10);
fprintf(fid, '  ]\n');

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
