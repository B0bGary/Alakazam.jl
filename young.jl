mutable struct YoungTableau
    row_lengths::Vector{Integer}
    index_assignment::Vector{Integer}
    tableau::Matrix{Integer}
    adj_matrix::Matrix{Integer}
end
    function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64})
        YoungTableau(row_lengths,index_assignment,maximum(index_assignment))
    end


    function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64},tot_inds::Int64)


        #create tabelau in matri form
        matrix =[-1 for _ in 1:size(row_lengths)[1], _ in 1:row_lengths[1]]
        cum=1
        for row in 1:size(row_lengths)[1]
            for i in 1:row_lengths[row]
                matrix[row,i]=index_assignment[cum]
                cum+=1

            end
        end
        println(matrix)
        println("here")


        #add symmetric perms
        adjm =[0 for _ in 1:tot_inds, _ in 1:tot_inds]
        init=-1
        for row in 1:size(matrix)[1]
            init*=-1
            cur=matrix[row,1]


            if cur!=-1
            for col in 1:size(matrix)[2]
                nex=matrix[row,col]
                if nex==-1
                    break
                else
                    last_col=col
                end



                while last_col>0
                    cr=matrix[row,last_col]
                    for cc in 1:last_col
                        nex=matrix[row,cc]
                        adjm[cr,nex]=1
                        adjm[nex,cr]=1
                    end
                    last_col-=1
                end
            end
            end
        

        for col in 1:size(matrix)[2]
            entry=matrix[row,col]
            if entry==-1
                break
            else
                for row2 in row+1:size(matrix)[1]
                    for col2 in 1:size(matrix)[2] 
                        nex=matrix[row2,col2]
                        if nex==-1
                            break
                        else
                            adjm[entry,nex]=-1
                            adjm[nex,entry]=-1
                        end
                    end
                end
            end
        end
    end


        println(matrix)
        println(adjm)





        YoungTableau(row_lengths, index_assignment,matrix,adjm)
    end