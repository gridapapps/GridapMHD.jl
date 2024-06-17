using NLsolve

function mesh_params(;Ha=100.0, #Hartmann number refered to the outlet channel
                     N_Ha=20,   #Total number of nodes along the Hartmann direction
                     N_side=20, #Total number of nodes along the side direction
                     n_Ha=4,    #Required number of nodes in the Hartmann BL
                     n_side=5,  #Required number of nodes in the side BL
		     Z=4.0, 	#Expansion ratio
                     β=0.2      #Aspect ratio of the outlet channel
                     )

#Equation for solving the geometric ratio nessary to have n cells in a boundary layer of ξ thickness

       f(r;a,ξ,N,n) = (1-r[1]^N)/(1-r[1]^n)-a/ξ
       f_Ha(r) = f(r; a = (Z-1)/(2*Z), ξ = 1/Ha, N = N_Ha/2, n = n_Ha)
       f_side(r) = f(r; a = β,ξ = 1/((Ha*Z)^0.5), N = N_side/2, n = n_side)

       sol_Ha = nlsolve(f_Ha,[2.0])
       sol_side = nlsolve(f_side,[2.0])

       params = [sol_Ha.zero[1],sol_side.zero[1]]
end

function R_calc(;   N=20,   #Total number of nodes along the segment
                    Dx=0.1, #Distance of the first node
		    L=1.0   #Total lenght of the segment
                     )

#Equation for solving the geometric ratio nessary to have the appropiate first thickness

       f(r;_N,_Dx) =_Dx*(1-r[1]^_N)/(1-r[1])-L
       fs(r) = f(r;_N=N,_Dx=Dx)
       sol = nlsolve(fs,[2.0])
       R= sol.zero[1]
       R
end
