// Gmsh project created on Tue Jan 16 14:55:13 2023

// (Externally written) Mesh parameters
Ha=1000.0;
N_Ha=22;
N_s=18;
r=4.0;
betta=0.2;
L_in=1.0;
L_out=2.0;
N_1=50;
N_2=80;
R_1=1.0760308637925087;
R_2=1.0251281583778045;
dx=0.05;


// Points to generate the channel source edge (y = -1 and z = -betta)

//Outlet channel 
Point(1) = {0, -1, -betta, 1.0};
Point(2) = {L_out, -1, -betta, 1.0};

Line(1) = {1,2};

//Inlet channel (for the moment the expansion takes place only in the field direction)

Point(3) = {-L_in, -1/r, -betta, 1.0};
Point(4) = {0, -1/r, -betta, 1.0};

Line(2) = {3,4};

//Uniform mesh along the flow direction

Transfinite Curve {1} = N_1 Using Progression R_1; 
Transfinite Curve {-2} = N_2 Using Progression R_2;


//Start saving the physical points in list for the tags
wall_lines[] = {1,2};
wall_points[] = {1, 2, 3, 4};

//------------------------------------------------------------------
//Extrusion in the direction perpendicular to B according to Smolentsev formula

one[0] = 1;
uniform[0] = 1/N_s;
a = (Ha^0.5/(Ha^0.5-1))^0.5;
c = (a + 1)/(a - 1);

s[0] = c^(2*(uniform[0]-0.5)); 
layer[0] = ((a+1)*s[0] - a + 1)/(2*(1+s[0]));

For i In {1:N_s-1}
   one[i] = 1;
   uniform[i] = uniform[i-1]+1/N_s;
   s[i] = c^(2*(uniform[i]-0.5)); 
   layer[i] = ((a+1)*s[i] - a + 1)/(2*(1+s[i]));
EndFor

extruded_entities_side[] = Extrude {0.0,0.0,2*betta} {Line{1,2};Layers{one[],layer[]}; Recombine;};
//-----------------------------------------------------------------------------

//Physical tags from the extruded list
wall_surfaces[] = {extruded_entities_side[1],extruded_entities_side[5]};

wall_lines[2] = extruded_entities_side[0];
wall_lines[3] = extruded_entities_side[2];
wall_lines[4] = -extruded_entities_side[3];
wall_lines[5] = extruded_entities_side[4];
wall_lines[6] = extruded_entities_side[6];
wall_lines[7] = -extruded_entities_side[7];

// From the GUI
For i In {5:8}
wall_points[i-1] = i;
EndFor

//----------------------------------------------------------------
//Frist extrusion along the direction of the B field (right part);

one[0] = 1;
uniform[0] = 1/N_Ha;
a = (Ha/(Ha-1))^0.5;
c = (a + 1)/(a - 1);

s[0] = c^(2*(uniform[0]-0.5)); 
layer[0] = ((a+1)*s[0] - a + 1)/(2*(1+s[0]));

For i In {1:N_Ha-1}
   one[i] = 1;
   uniform[i] = uniform[i-1]+1/N_Ha;
   s[i] = c^(2*(uniform[i]-0.5)); 
   layer[i] = ((a+1)*s[i] - a + 1)/(2*(1+s[i]));
EndFor


extruded_entities_Ha[] = Extrude {0.0,(r-1)/r,0.0} {Surface{extruded_entities_side[1]};Layers{one[],layer[]}; Recombine;};

//----------------------------------------------------------------

//Physical tags from the extruded list
volumes[] = extruded_entities_Ha[1]; 

wall_surfaces[2] = extruded_entities_Ha[2]; 
wall_surfaces[3] = extruded_entities_Ha[4]; 
wall_surfaces[4] = extruded_entities_Ha[5]; 

outlet_surfaces[0] = extruded_entities_Ha[3]; 

// From the GUI
wall_lines_2[] = {12, 14, 17, 18, 22, 26};
wall_points_2[] = {10, 14};

outlet_lines[0] = 13;

//-----------------------------------------------------------------------------
//Second extrusion along the direction of the B field (central part and inlet channel);

extruded_entities_Ha[] = Extrude {0.0,2/r,0.0} {Surface{extruded_entities_Ha[0],extruded_entities_side[5]};Layers{one[],layer[]}; Recombine;};
//--------------------------------------------------------------

//Physical tags from the extruded list
volumes[1] = extruded_entities_Ha[1];
volumes[2] = extruded_entities_Ha[7];

wall_surfaces[6] = extruded_entities_Ha[2];
wall_surfaces[7] = extruded_entities_Ha[4];
outlet_surfaces[1] = extruded_entities_Ha[3];

wall_surfaces[8] = extruded_entities_Ha[6];
wall_surfaces[9] = extruded_entities_Ha[8];
wall_surfaces[10] = extruded_entities_Ha[10];
inlet_surfaces[0] = extruded_entities_Ha[11];

// From the GUI
wall_lines_3[] = {34, 36, 37, 39, 40, 44, 48, 56, 58, 59, 61, 70};
wall_points_3[] = {15, 16, 20, 25, 24, 34};

outlet_lines[1] = 35;

//-----------------------------------------------------------------------------
//Third extrusion along the direction of the B field (left part and inlet channel);

extruded_entities_Ha[] = Extrude {0.0,(r-1)/r,0.0} {Surface{extruded_entities_Ha[0]};Layers{one[],layer[]}; Recombine;};
//-----------------------------------------------------------------------------

//Physical tags from the extruded list
volumes[3] = extruded_entities_Ha[1];

wall_surfaces[11] = extruded_entities_Ha[0];
wall_surfaces[12] = extruded_entities_Ha[2];
wall_surfaces[13] = extruded_entities_Ha[4];
wall_surfaces[14] = extruded_entities_Ha[5];

outlet_surfaces[2] = extruded_entities_Ha[3];

wall_lines_4[] = {78, 79, 80, 81, 83, 84, 88, 92};
wall_points_4[] = {35, 36, 40, 44};

//-------------Definition of physical tags from the lists---------------

Physical Point ("wall") = {wall_points[], wall_points_2[], wall_points_3[], wall_points_4[]};
Physical Line ("wall") = {wall_lines[], wall_lines_2[], wall_lines_3[], wall_lines_4[]};
Physical Surface ("wall") = {wall_surfaces[]};

Physical Line ("outlet") = {outlet_lines[]};
Physical Surface ("outlet") = {outlet_surfaces[]};

Physical Surface ("inlet") = {inlet_surfaces[]};

Physical Volume("PbLi") = {volumes[]};
