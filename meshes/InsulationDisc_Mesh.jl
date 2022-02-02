using Gridap
using GridapMHD

# Rectangular insulated channel with a discontinuity in the insulation

#----------Problem setting------- 

#Dimensions of the channel ((L_in+L_out)*2a*2b) 2tGap is the thickness of the insulation discontinuity 
  a = 1.0     
  b = 0.625*a
  tGap = 0.25*a
  LIn = 10*a
  LOut = 20*a
  tGap_ext = a             #Thickness for the concentration of nodes next to the gap 

#Mesh parameters (hexas)
  nx = 300
  ny = 40
  nz = 25
  
  nIn  = 80		#Node in the inlet channel (even number)
  nGap = 60   		#Node concentration in the expected 3D region(even number)
  strech_c = 1.2 	#Streching towards the gap

# The Hartmann number determines the streching towards the boundary layers according with the work of Smolentsev

  Ha = 1000

#-----------------Functions for mesh construction------------------  
 
function strechMHD(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),factor=(1.0,1.0,1.0),
  dirs=(1,2,3))
ncoord = collect(coord.data)
for (i,dir) in enumerate(dirs)
  ξ0 = domain[i*2-1]
  ξ1 = domain[i*2]
  l =  ξ1 - ξ0
  c = (factor[i] + 1)/(factor[i] - 1)

  if l > 0
        if ξ0 <= coord[dir] <= ξ1
           ξx = (coord[dir] - ξ0)/l                     #ξx from 0 to 1 uniformly distributed
           ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)    #ξx streched from 0 to 1 towards 1
           ncoord[dir] =  ξx_streched*l + ξ0            #coords streched towards ξ1
        end
  else
        if ξ1 <= coord[dir] <= ξ0
           ξx = (coord[dir] - ξ0)/l                        #ξx from 0 to 1 uniformly distributed
           ξx_streched = factor[i]*(c^ξx-1)/(1+c^ξx)       #ξx streched from 0 to 1 towards 1
           ncoord[dir] =  ξx_streched*l + ξ0               #coords streched towards ξ1

        end
  end
end
VectorValue(ncoord)
end

function ChangeDensity(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),subDomain=(0.0,1.0,0.0,1.0,0.0,1.0), 
			     nodesTot=(1.0,1.0,1.0), nodesSub=(1.0,1.0,1.0), dirs=(1,2,3))
ncoord = collect(coord.data)
for (i,dir) in enumerate(dirs)
  ξ0 = domain[i*2-1]
  ξ1 = domain[i*2]
  Ltot =  ξ1 - ξ0
  Lsub = subDomain[i*2] - subDomain[i*2-1]
  
  alpha = (Lsub/Ltot)*(nodesTot[i]/nodesSub[i])
  betta = ((Ltot-Lsub)/Ltot)*(nodesTot[i]/(nodesTot[i]-nodesSub[i]))

  if Ltot != Lsub

     if ξ0 <= coord[dir] <= ξ1
       
        if coord[dir] <= (Lsub/alpha + ξ0)
	   
           ncoord[dir] = alpha*coord[dir] + ξ0*(1-alpha)
	
	else
	   
   	   ncoord[dir] = betta*coord[dir] + ξ0*(1-betta) + Lsub*(1-betta/alpha)

        end
     else
	  ncoord[dir] = coord[dir]

     end
  end
end
VectorValue(ncoord)
end

#------------Mesh and geometry-------------------

  # Background uniform mesh
  domain = (-LIn/a,LOut/a,-1.0,1.0,-b/a,b/a)
  partition = (nx,ny,nz)  


  #Streching and node concentration
  strech_a = sqrt(Ha/(Ha-1))
  strech_b = sqrt(sqrt(Ha)/(sqrt(Ha)-1))

function map1(coord)
   ncoord = ChangeDensity(coord; domain=(-LIn/a, LOut/a), subDomain=(-LIn/a, -tGap_ext/a), nodesTot=(nx,), nodesSub=(nIn,), dirs=(1,))
   ncoord = ChangeDensity(ncoord;domain=(-tGap_ext/a,LOut/a),subDomain=(-tGap_ext/a,tGap_ext/a), nodesTot=((nx-nIn),), nodesSub=(nGap,), dirs=(1,))
   ncoord = strechMHD(ncoord, domain=(-LIn/a,-tGap_ext/a,0.0,-a,0.0,-b/a),factor=(strech_c,strech_a,strech_b),dirs=(1,2,3))
   ncoord = strechMHD(ncoord, domain=( LOut/a, tGap_ext/a,0.0, a,0.0, b/a),factor=(strech_c,strech_a,strech_b),dirs=(1,2,3))
   ncoord
end

  model = CartesianDiscreteModel(domain,partition,map=map1)

#-----------Boundary Conditions--------------------------------

# Body force
  f_u = VectorValue(0.0,0.0,0.0) #Source term

  # Parabolic or constant profile at the inlet
function g_u(coord)
  x = coord[1]
  y = coord[2]
  z = coord[3]
  tol = 1e-6
  if x < (-LIn/a + tol)
   #  return VectorValue((9/4)*(a/b)^2 * (b/a-z) * (b/a+z) * (1-y) * (1+y), 0.0, 0.0)
     return VectorValue(1.0, 0.0, 0.0)
  else
    return VectorValue(0.0, 0.0, 0.0)
  end
end 


 g_j = VectorValue(0.0,0.0,0.0)


# Filter to select all walls
function walls(coord)
  x = coord[1]
  y = coord[2]
  z = coord[3]
  tol = 1e-6
  if (x < -tGap/a - tol) | (x > tGap/a + tol)
     if (z < -b/a + tol) | (z > b/a - tol)
    	return true
     elseif (y < -1.0 + tol ) | ( y > 1.0 - tol)
    	return true
     end
  end
  false
end

# Filter to select the wall gap
function gap(coord)
  x = coord[1]
  y = coord[2]
  z = coord[3]
  tol = 1e-6
  if (x > -tGap/a - tol) & (x < tGap/a + tol)
     if (z < -b/a + tol) | (z > b/a - tol)
        return true
     elseif (y < -1.0 + tol ) | ( y > 1.0 - tol)
        return true
     end
  end
  false
end

# Filter to select inflow boundary
function inflow(coord)
  x = coord[1]
  y = coord[2]
  z = coord[3]
  tol = 1e-6
  if (x < -LIn/a + tol)
    return true
  end
  false
end

# Filter to select outflow boundary
function outflow(coord)
  x = coord[1]
  y = coord[2]
  z = coord[3]
  tol = 1e-6
  if (x > LOut/a - tol)
    return true
  end
  false
end

# Name inlet, outlet, and walls element sets.
add_entity!(model,walls,"walls")
add_entity!(model,gap,"gap")
add_entity!(model,inflow,"inlet")
add_entity!(model,outflow,"outlet")

writevtk(model,"Mesh")


