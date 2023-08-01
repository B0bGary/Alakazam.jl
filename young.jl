mutable struct YoungTableau
    row_lengths::Vector{Integer}
    index_assignment::Vector{Integer}
    tableau::Matrix{Integer}

    function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64})

        matrix =[-1 for _ in 1:size(row_lengths)[1], _ in 1:row_lengths[1]]
        cum=1
        for row in 1:size(row_lengths)[1]
            for i in 1:row_lengths[row]
                matrix[row,i]=index_assignment[cum]
                cum+=1

            end
        end

        new(row_lengths, index_assignment,matrix)
    end
end

