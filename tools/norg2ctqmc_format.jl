using DelimitedFiles
using Printf

# 假设的文件处理函数
function process_file(filepath, old_prefix, new_prefix)
    # 读取文件
    matrix_data = readdlm(filepath)
    matrix_data = matrix_data[2:8194, :]
    new_filename = replace(filepath, old_prefix => new_prefix)

    # 初始化目标数组
    result = zeros(5, 8193, 4)  
    band = 5

    # 将数据填充到目标数组
    for i in 1:8193
        # 第一组
        result[1, i, 1] = 1
        result[1, i, 2:3] = matrix_data[i, 1:2]
        result[1, i, 4] = matrix_data[i, 28]
        
        # 第二组
        result[2, i, 1] = 2
        result[2, i, 2] = matrix_data[i, 1]
        result[2, i, 3] = matrix_data[i, 8]
        result[2, i, 4] = matrix_data[i, 34]
        
        # 第三组
        result[3, i, 1] = 3
        result[3, i, 2] = matrix_data[i, 1]
        result[3, i, 3] = matrix_data[i, 14]
        result[3, i, 4] = matrix_data[i, 40]
        
        # 第四组
        result[4, i, 1] = 4
        result[4, i, 2] = matrix_data[i, 1]
        result[4, i, 3] = matrix_data[i, 20]
        result[4, i, 4] = matrix_data[i, 46]
        
        # 第五组
        result[5, i, 1] = 5
        result[5, i, 2] = matrix_data[i, 1]
        result[5, i, 3] = matrix_data[i, 26]
        result[5, i, 4] = matrix_data[i, 52]
    end

    # # 返回结果
    # result
    # 打开文件以追加模式
    # 打开文件以追加模式
    open(new_filename, "w") do f
        for group in 1:band
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

    open(new_filename, "a") do f
        for group in 1:band
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
end

function process_and_rename_files()
    for file in readdir(".")
        if occursin("gfimp.txt.", file)
            # 处理文件
            process_file(file, "gfimp.txt.", "solver.grn.dat.")
        end
    end
end

# 调用函数
process_and_rename_files()





