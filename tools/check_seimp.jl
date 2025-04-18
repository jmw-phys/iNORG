using DelimitedFiles

norb = 3

# 获取当前目录下的所有文件名
all_files = readdir()

# 筛选出所有以 "seimp.txt." 开头的文件
target_files = filter(filename -> startswith(filename, "seimp.txt"), all_files)

# 对筛选出的文件进行操作
for filename in target_files
    m = readdlm(filename)
    println("Processing file: ", filename)
    for j = (norb*norb) + 3 : (norb*norb) * 2 + 2
        for i = 2 : 10
            if m[i, j] > 0 
                println(filename, "  ", m[i, (norb*norb) + 2], "  ", m[1, j], "  ", m[i, j])  
            end
        end
    end
end
