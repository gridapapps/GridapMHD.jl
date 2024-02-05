using NLsolve

function mesh_params(;Ha=100.0, #Hartmann number
                     N_Ha=20,   #Total number of nodes along the Hartmann direction
                     N_side=20, #Total number of nodes along the Hartmann direction
                     n_Ha=4,    #Required number of nodes in the Hartmann BL
                     n_side=5,  #Required number of nodes in the side BL 
                     β=0.25     #Aspect ratio of the inlet channel (expansion case)
                     )
                     
#Equation for solving the geometric ratio nessary to have n cells in a boundary layer of ξ thickness       

       f(r;ξ,N,n) = r[1]^n*(1-ξ*r[1]^(N-n))-ξ-1.0 
       
       f_Ha(r) = f(r; ξ = 1/Ha, N = N_Ha, n = n_Ha)
       f_side(r) = f(r; ξ = β/(Ha^0.5), N = N_side, n = n_side)  
       
       sol_Ha = nlsolve(f_Ha,[2.0])
       sol_side = nlsolve(f_side,[2.0])
                    
       params = [sol_Ha.zero[1],sol_side.zero[1]]
end
