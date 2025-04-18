# 假设 matrix_data 是您的 9549x7 的数据矩阵
# matrix_data = ...
using DelimitedFiles
using Printf
matrix_data = readdlm("Delta.inp")

# 初始化目标数组
result = zeros(3, 8193, 4)  # 注意这里的维度更改为 4

# 将数据填充到目标数组
for i in 1:8193
    # 第一组
    result[1, i, 1] = 1
    result[1, i, 2:4] = matrix_data[i, 1:3]
    
    # 第二组
    result[2, i, 1] = 2
    result[2, i, 2] = matrix_data[i, 1]
    # result[2, i, 3:4] = matrix_data[i, 4:5]
    result[2, i, 3:4] = matrix_data[i, 2:3]
    
    # 第三组
    result[3, i, 1] = 3
    result[3, i, 2] = matrix_data[i, 1]
    # result[3, i, 3:4] = matrix_data[i, 6:7]
    result[3, i, 3:4] = matrix_data[i, 2:3]
end

# # 返回结果
# result
# 打开文件以追加模式
# 打开文件以追加模式
open("solver.hyb.in", "w") do f
    for group in 1:3
        for i in 1:8193
            # 将数组中的每一行格式化为字符串并写入文件
            line = @sprintf("     %d      %.10f      %.10f      %.10f      0.00000000      0.00000000\n", 
                            Int(result[group, i, 1]), result[group, i, 2], result[group, i, 3], result[group, i, 4])
            write(f, line)
        end
        # 每组数据之间添加一个空行以区分
        write(f, "\n\n")
    end
end

open("solver.hyb.in", "a") do f
    for group in 1:3
        for i in 1:8193
            # 将数组中的每一行格式化为字符串并写入文件
            line = @sprintf("     %d      %.8f      %.8f      %.8f      0.00000000      0.00000000\n", 
                            Int(result[group, i, 1]+3), result[group, i, 2], result[group, i, 3], result[group, i, 4])
            write(f, line)
        end
        # 每组数据之间添加一个空行以区分
        write(f, "\n\n")
    end
end