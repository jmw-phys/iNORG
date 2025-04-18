
using DelimitedFiles
using Plots

# cutoff_value = parse(Int, match(r"--cutoff=(\d+)", ARGS[1])[1])
# x_value = parse(Float64, match(r"--x=([0-9\.]+)", ARGS[1])[1])
# println(x_value)
b = readdir()
# xsize = 59999
xsize = 1
accurt = 0.0000001

for a in b
    if (a != "code" && a != "output" ) 
        m = readdlm(a)
        println("test the ",a)
        no = reverse(m[1:size(m)[1],4])
        no = no[1:size(no,1)-1]
        mstate = zeros(xsize,2)
        left = 0
        for cutoff=1:xsize
            sum = 0.0
            for j in eachindex(no)
                sum += no[j]
                # cutoff = 0.5
                # cutoff = x_value
                left = size(no,1) - j
                if(sum > cutoff*accurt) 
                    println("Cut for if sum = ", cutoff*accurt, " of whole state.")
                    println(j," The Full size ",size(no,1)," left size ",left," left ",100*(left)/size(no,1),"% ")
                    println("first and last element",no[1],"  ",no[size(no,1)])
                    mstate[cutoff,1] = cutoff*accurt
                    mstate[cutoff,2] = left
                    break
                end
            end
            if(left < 1e3) 
                # xsize = left
                mstate = mstate[1:cutoff,:]
                break
            end
        end
        println(mstate)
        mstate[:, 1] = round.(mstate[:, 1], digits=6)  # 保留第一列小数点后两位
        if !isdir("output")
            mkdir("output")
        end        
        dim = size(no,1)
        writedlm("output/$(dim)output_$(accurt)left$(left).txt", mstate, repeat(' ', 6))  # 用制表符分隔数据
        # writedlm("output.txt", round.(mstate, digits=2), '\t', false, "%0.2f\t%d")
    end
    # # 生成一个示例矩阵
    # mstate = hcat(collect(1:xsize)', rand(xsize))
    
    # # 将第二列数据取对数
    # mstate[:, 2] = log10.(mstate[:, 2])
    
    # # 绘制散点图，使用对数坐标轴
    # scatter(mstate[:, 1]', mstate[:, 2], yaxis=:log10)
    # break
end

