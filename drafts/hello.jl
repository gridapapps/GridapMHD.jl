using PartitionedArrays
using GridapMHD

np=4
with_mpi() do distribute
    ranks = distribute(LinearIndices((np,)))
    map(ranks) do rank
       println("I am proc $rank of $np")
    end  
end
