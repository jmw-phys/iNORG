using DelimitedFiles
using Printf
using Statistics  # 引入统计库用于计算标准差

# 修改函数以读取目录中所有文件并计算矩阵的算术平均和跨文件的标准差
function read_ctqmc_data(path)
    files = readdir(path)
    all_data = []
    col3_data = []  # 存储所有文件的第3列数据
    col4_data = []  # 存储所有文件的第4列数据

    for file in files
        filepath = joinpath(path, file)
        data = readdlm(filepath, Float64)
        push!(all_data, data)

        # 收集第3列和第4列数据
        push!(col3_data, data[:, 3])
        push!(col4_data, data[:, 4])
    end

    if length(all_data) == 0
        return nothing
    end

    mean_matrix = zeros(size(all_data[1]))
    for matrix in all_data
        mean_matrix .+= matrix
    end
    mean_matrix ./= length(all_data)

    # 计算跨文件的第3列和第4列标准差
    std_devs_col3 = std(hcat(col3_data...), dims=2)
    std_devs_col4 = std(hcat(col4_data...), dims=2)

    # 打印标准差以供检查
    println("Standard deviation for col3: ", std_devs_col3[1,:])
    println("Standard deviation for col4: ", std_devs_col4[1,:])

    # 将标准差放入最终矩阵的第4列和第5列
    mean_matrix[:, 5] = std_devs_col3
    mean_matrix[:, 6] = std_devs_col4
    
    return mean_matrix
end

# Function to group data based on the first column
function group_data(data, nband)
    grouped_data = Dict()
    for row in eachrow(data)
        group_key = row[1]
        values = convert(Vector{Float64}, row[2:end])
        if group_key > nband break end
        push!(get!(grouped_data, group_key, []), values)
    end

    # Convert to Matrix{Float64}
    for (key, value) in grouped_data
        grouped_data[key] = hcat(value...)'
    end
    
    return grouped_data
end

# Function to save data to .data files
function save_data(grouped_data, x)
    for (key, value_matrix) in grouped_data
        # Save only the first x rows
        partial_matrix = value_matrix[1:min(x, end), :]
        filename = "band_$(Int(key)).data"
        open(filename, "w") do io
            for row in eachrow(partial_matrix)
                formatted_row = join([@sprintf("%18.8f", value) for value in row], " ")
                write(io, formatted_row, "\n")
            end
        end
    end
end


# Run the main function
filepath = "solver.grn"  # Replace with your file path
lent = 102  # Save only the first x rows for each group
nband = 5
# main(filepath, lent)
data = read_ctqmc_data(filepath)
grouped_data = group_data(data, nband)
save_data(grouped_data, lent)