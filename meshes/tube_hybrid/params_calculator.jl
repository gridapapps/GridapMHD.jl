using NLsolve

function mesh_params(;Ha=100.0, #Hartmann number
                     N = 20,    #Total number of nodes along the Hartmann direction
                     n = 4,     #Required number of nodes in the Hartmann BL
                     a = 0.5,   #Total lenght of the nodalized curve 
                     )
                     
#Equation for solving the geometric ratio nessary to have n cells in a boundary layer of ξ thickness in a curve with lenght a      

       f(r;ξ,N,n,a) = (1-r[1]^N)/(1-r[1]^n)-a/ξ 
       
       f_Ha(r) = f(r; ξ = 1/Ha, N = N, n = n, a = a)
       
       sol = nlsolve(f_Ha,[2.0])
                           
       sol.zero[1]
end
