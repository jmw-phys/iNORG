#!/usr/bin/env julia

push!(LOAD_PATH, ENV["ACFLOW_HOME"])
using ACFlow

using DelimitedFiles
using Printf

# 获取当前目录下的所有文件
files = readdir()

# 定义需要提取的列组合和对应的文件名前缀
# columns_and_prefixes = [([1, 2, 7], "1B_"), ([1, 5, 10], "2B_")] #two orbitals
columns_and_prefixes = [([1, 2, 12], "1B_"), ([1, 6, 16], "2B_"), ([1, 10, 20], "3B_")] #three orbitals
# columns_and_prefixes = [([1, 2, 5], "WB_")]

# 确保 "DOS" 文件夹存在
dos_dir = "DOS"
if !isdir(dos_dir)
    mkdir(dos_dir)
end

for file in files
    # 检查文件是否是 .txt 文件
    if occursin(".txt", file)
        # 读取文件
        data = readdlm(file)

        for (columns, prefix) in columns_and_prefixes
            # 提取指定的列
            new_data = data[2:61, columns]

            # 添加一列全为0.0001的数据
            rows, cols = size(new_data)
            new_column = fill(0.0001, rows)
            new_data = hcat(new_data, new_column)

            # 格式化数据，保留6位小数
            formatted_data = [@sprintf("%12.6f", val) for val in new_data]

            # 将结果写入新的文件
            # 替换文件后缀为 ".data"
            new_file_name = replace(file, ".txt" => ".data")
            new_file_name = prefix * new_file_name
            # 放入 "DOS" 文件夹中
            new_file_name = joinpath(dos_dir, new_file_name)
            writedlm(new_file_name, formatted_data)

            # MaxEnt solver
            welcome()

            # Setup parameters
            B = Dict{String,Any}(
                "finput" => new_file_name,
                "mtype"  => "flat",
                "mesh"   => "linear",
                "ngrid"  => 60,
                "nmesh"  => 2001,
                "wmax"   => 5.0,
                "wmin"   => -5.0,
                "beta"   => 314.159265358979323,
            )

            S = Dict{String,Any}(
                "nalph"  => 15,
                "alpha"  => 1e15,
            )

            setup_param(B, S)

            # Call the solver
            mesh, Aout, Gout = solve(read_data())

            # Move and rename "Aout.data"
            old_Aout_path = "Aout.data"
            # new_Aout_path = joinpath("DOS", prefix * "Aout.data")
            new_Aout_path = joinpath(dos_dir, prefix * replace(file, "mb.gfimp.txt" => "Aout.txt"))
            mv(old_Aout_path, new_Aout_path, force=true) # force=true 会覆盖同名的目标文件

            # Remove other generated files
            for generated_file in ["chi2.data", "Gout.data", "model.data", "repr.data"]
                rm(generated_file, force=true) # force=true 会忽略不存在的文件
            end
        end
    end
end
