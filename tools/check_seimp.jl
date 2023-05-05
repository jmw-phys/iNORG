
using DelimitedFiles

# println("The out put is in: "*pwd()*"/out")# for customers.

# mkdir("out")
norb = 2
b = readdir()
for a in b
    if (a != "code") 
        m = readdlm(a)
        println("test the ",a)
        for j=(norb*norb)+3:(norb*norb)*2+2
            for i=2:10
                if(m[i, j]>0) println(a,"  ",m[i, 27],"  ",m[1, j],"  ",m[i, j])  end
            end
        end
    end
end