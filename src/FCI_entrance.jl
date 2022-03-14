# Change of the electrical insulation at the entrance of an FCI
# Validation problem with respect to the experiment published by L. Bühler et al. "Experimental investigation of liquid metal MHD flow entering a flow channel insert"
# Fusion Engineering and Design 154(2020):111484

using Gridap
using GridapGmsh
using GridapDistributed
using PartitionedArrays
using GridapMHD

#----------Problem setting------- 


function FCI_entrance(;
   Ha = 1.0,	  			#Hartmann number
   Re = 1, 	  			#Reynolds number
   c_w = 0.0727,  			#Wall conductance ratio 
   c_FCI = 0.00476,		 	#FCI conductance ratio
   B_dir = VectorValue(0.0,1.0,0.0), 	#Direction of the B-field (norm =1)  
   α_w = 100,	  			#Penalty for the wall boundary condition
   α_FCI =100, 	  			#Penalty for the FCI boundary condition
   R = 1,				#Radious of the pipe
   L = 6,	 			#Half length of the channel (normalized by R)	   
   ρ = 1, 		      		#Density (Kg/m3)
   ν = 1,      				#Kinematic viscosity (Pa·s)
   σ = 1,   				#Electrical conductivity (S/m)
   np = 1				#Number of parts for parallel computations 
)

# Derived variables 
  
  N = Ha^2/Re			   #Interaction parameter or Stuart number
  
  U0 = Re*ν/R			   #Velocity scale	
  B0 = Ha*sqrt((ρ*ν)/σ)/R 	   #Module of the external B-field	

#------------Mesh and geometry------------------------

  #This is for parallel id
  parts = get_part_ids(sequential,1)

  model = GmshDiscreteModel("./geometries/FCI_Entrance.msh") 
  Ω = Triangulation(model)

  t = PTimer(parts,verbose=true)

#-----------Boundary Conditions-----------------------

# Body force
  
  z = VectorValue(0.0, 0.0, 0.0)

  # Parabolic profile at the inlet
  function g_u(coord,R,L)
     r = sqrt(coord[1]^2 + coord[2]^2)
     zz = coord[3]
     tol = 1e-6
     if zz < (-L/R + tol)
         return VectorValue(0.0, 0.0, U0*(R^2-r^2))
     else
    	 return VectorValue(0.0, 0.0, 0.0)
     end
  end
  
  #Thin wall boundaries

  Γ_w   = Boundary(model, tags="wall")
  Γ_FCI = Boundary(model, tags="FCI")

#-------Dictionary of parameters for the solver-----------

 params = Dict(
    :ptimer=>t,
    :debug=>true,
    :fluid=>Dict(
      :domain=>model,
      :α=>1,
      :β=>1/Re,
      :γ=>N,
      :u=>Dict(
        :tags=>["wall","FCI","inlet"],
        :values=>[z,z,g_u]),
      :j=>Dict(
        :tags=>["inlet", "outlet"],
        :values=>[z,z]),
#      :thin_wall=>[Dict(:domain=>Γ_w, 	 :cw=>c_w,   :τ=>τ_w),
#		    Dict(:domain=>Γ_FCI, :cw=>c_FCI, :τ=>τ_FCI)], 
      :f=>z,
      :B=>B_dir, 
    ),
  )

#-----------------Solver call-------------------------

xh = main(params)

uh, ph, jh, φh = xh
divj = ∇⋅jh

writevtk(Ω,"Hydrodynamic",cellfields=["u"=>uh,"p"=>ph,"j"=>jh,"phi"=>φh,"divj"=>divj])


end
