%% 颜色生成函数
function c = get_color(idx)
    colors = lines(7);
    c = colors(mod(idx-1,7)+1,:);
end