function result = my_prod(arr)
    % 输入：arr 是一个数组
    % 输出：result 是数组中所有元素的乘积

    result = 1; % 初始化为 1
    for i = 1:length(arr)
        result = result * arr(i);
    end
end