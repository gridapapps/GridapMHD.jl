// Gmsh project created on Tue Jan 16 14:55:13 2023

// (Externally written) Mesh parameters
Ha=1000.0;
N_Ha=6;
N_s=6;
exp_Ha=0.6;
exp_s=0.5;
r=4.0;
betta=1.0;
L_in=8.0;
L_out=8.0;
t=0.028;
N_1=10;
N_2=10;
N_w=1;
R_1=1.4278802278903293;
R_2=1.4278802278903293;


// Points to generate the channel source edge

//Outlet channel  (y = -1-t and z = -betta-t)
Point(1) = {-t, -1-t, -betta-t, 1.0};
Point(2) = {0, -1-t, -betta-t, 1.0};
Point(3) = {L_out, -1-t, -betta-t, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};

//Inlet channel (y = -1/r-t, z = -betta-t)

Point(4) = {-L_in, -1/r-t, -betta-t, 1.0};
Point(5) = {-t, -1/r-t, -betta-t, 1.0};
//Point(6) = {0, -1/r-t, -betta-t, 1.0};

Line(3) = {4,5};
//Line(4) = {5,6};

//Start saving the physical points in list for the tags
ext_lines[] = {1, 2, 3};
ext_points[] = {1, 2, 3, 4, 5};

//-----------------------------------------------------
//Mesh along the flow direction (refinement towards the expansion)

Transfinite Curve {-3} = N_1 Using Progression R_1; 
Transfinite Curve {2} = N_2 Using Progression R_2;

//-------------------------------------------------------------
//Uniform wall mesh in the x direction

Transfinite Curve {1} = N_w; 

//---------------------------------------------------------
// Uniform extrusion in the direction perpendicular to B (wall)

extruded_entities_ext[] = Extrude {0.0,0.0,t} {Line{1,2,3};Layers{N_w}; Recombine;};

//For i In {0:11}
//	Printf("%f",extruded_entities_ext[i]);
//EndFor

//Physical tags from the extruded list
ext_surfaces[] = {extruded_entities_ext[1],extruded_entities_ext[5],extruded_entities_ext[9]};

ext_lines_2[0] = extruded_entities_ext[0];
k = 1;
For n In {1:2}
	For i In {4*n-2:4*n}
		k = k+1;
		ext_lines_2[k] = extruded_entities_ext[i];
//		Printf("%f",k);
	EndFor
EndFor
ext_lines_2[8] =  extruded_entities_ext[10];
ext_lines_2[9] =  extruded_entities_ext[11];

ext_lines[] ={ext_lines[],ext_lines_2[]};

//From GUI
ext_points[] = {ext_points[],6, 7,9,10,11};

//------------------------------------------------------------------
//Extrusion in the direction perpendicular to B according to Smolentsev formula (fluid)

one[0] = 1;
uniform[0] = 1/N_s;
a = (Ha^exp_s/(Ha^exp_s-1))^0.5;
c = (a + 1)/(a - 1);

s[0] = c^(2*(uniform[0]-0.5)); 
layer[0] = ((a+1)*s[0] - a + 1)/(2*(1+s[0]));

For i In {1:N_s-1}
   one[i] = 1;
   uniform[i] = uniform[i-1]+1/N_s;
   s[i] = c^(2*(uniform[i]-0.5)); 
   layer[i] = ((a+1)*s[i] - a + 1)/(2*(1+s[i]));
EndFor
//-----------------------------------------------------------------------------------

extruded_entities_side[] = Extrude {0.0,0.0,2*betta} {Line{4,8,12};Layers{one[],layer[]}; Recombine;};

//For i In {0:11}
//	Printf("%f",extruded_entities_side[i]);
//EndFor

//Physical tags from the extruded list
ext_surfaces[] = {ext_surfaces[],
			extruded_entities_side[1],extruded_entities_side[5],extruded_entities_side[9]};

ext_lines_3[0] = extruded_entities_side[0];
k = 1;
For n In {1:2}
	For i In {4*n-2:4*n}
		k = k+1;
		ext_lines_3[k] = extruded_entities_side[i];
	EndFor
EndFor
ext_lines_3[8] =  extruded_entities_side[10];
ext_lines_3[9] =  extruded_entities_side[11];

ext_lines[] ={ext_lines[],ext_lines_3[]};

// From the GUI
ext_points[] = {ext_points[],12,13,15,16,17};

//---------------------------------------------------------
// Uniform extrusion in the direction perpendicular to B (wall)

extruded_entities_ext[] = Extrude {0.0,0.0,t} {Line{16,20,24};Layers{N_w}; Recombine;};

//For i In {0:11}
//	Printf("%f",extruded_entities_ext[i]);
//EndFor

//Physical tags from the extruded list
ext_surfaces[] = {ext_surfaces[],
			extruded_entities_ext[1],extruded_entities_ext[5],extruded_entities_ext[9]};

ext_lines_4[0] = extruded_entities_ext[0];
k = 1;
For n In {1:2}
	For i In {4*n-2:4*n}
		k = k+1;
		ext_lines_4[k] = extruded_entities_ext[i];
//		Printf("%f",k);
	EndFor
EndFor
ext_lines_4[8] =  extruded_entities_ext[10];
ext_lines_4[9] =  extruded_entities_ext[11];

ext_lines[] ={ext_lines[],ext_lines_4[]};

//From GUI
ext_points[] = {ext_points[],18,19,21,22,23};
//---------------------------------------------------

//Frist extrusion along the direction of the B field (wall, outlet channel);

extruded_entities_wall[] = Extrude {0.0,t,0.0} {Surface{7,11,19,23,31,35};Layers{N_w}; Recombine;};

//For i In {0:35}
//	Printf("%f",extruded_entities_wall[i]);
//EndFor

//From list
wall_volumes[] = {extruded_entities_wall[1],extruded_entities_wall[7],extruded_entities_wall[13], extruded_entities_wall[19],extruded_entities_wall[25],extruded_entities_wall[31]};

ext_surfaces[] = {ext_surfaces[],
		extruded_entities_wall[2],extruded_entities_wall[8],extruded_entities_wall[5],extruded_entities_wall[17],extruded_entities_wall[29],extruded_entities_wall[28],extruded_entities_wall[34],extruded_entities_wall[9],extruded_entities_wall[21],extruded_entities_wall[33]};

//From GUI
int_surfaces[] = {127};
int_lines[] = {86,65,109,108}; 
int_points[] = {55,45,29,39};

ext_lines[] = {ext_lines[],88,131,132,153,63,161,139,69,143,99,47,46,55,44,41,152,117,73,64};
ext_points[] = {ext_points[],65,61,49,24,25,33,35,71};

//------------------------------------------------------------------
//Extrusion in the PbLi according to Smolentsev formula
//Second extrusion in the direction of the B-field

one[0] = 1;
uniform[0] = 1/N_Ha;
a = (Ha^exp_Ha/(Ha^exp_Ha-1))^0.5;
c = (a + 1)/(a - 1);

s[0] = c^(2*(uniform[0]-0.5)); 
layer[0] = ((a+1)*s[0] - a + 1)/(2*(1+s[0]));

For i In {1:N_Ha-1}
   one[i] = 1;
   uniform[i] = uniform[i-1]+1/N_Ha;
   s[i] = c^(2*(uniform[i]-0.5)); 
   layer[i] = ((a+1)*s[i] - a + 1)/(2*(1+s[i]));
EndFor
//---------------------------------------------------------------

extruded_entities_Ha[] = Extrude {0.0,(r-1)/r-t,0.0} {Surface{extruded_entities_wall[0],extruded_entities_wall[6],extruded_entities_wall[12],extruded_entities_wall[18],extruded_entities_wall[24],extruded_entities_wall[30]};
								  Layers{one[],layer[]}; Recombine;};

//For i In {0:35}
//	Printf("%f",extruded_entities_Ha[i]);
//EndFor

//Physical tags from the extruded list

wall_volumes[] = {wall_volumes[],
		extruded_entities_Ha[1],extruded_entities_Ha[7],extruded_entities_Ha[13],extruded_entities_Ha[25],extruded_entities_Ha[31]}; 
fluid_volumes[0] = extruded_entities_Ha[19];

int_surfaces[] = {int_surfaces[],
		extruded_entities_Ha[20], extruded_entities_Ha[22], extruded_entities_Ha[23]};
outlet_surfaces[0] = extruded_entities_Ha[21];
ext_surfaces[] = {ext_surfaces[],
		extruded_entities_Ha[2],extruded_entities_Ha[5],extruded_entities_Ha[8], extruded_entities_Ha[9],extruded_entities_Ha[17],extruded_entities_Ha[28], extruded_entities_Ha[29], extruded_entities_Ha[33],extruded_entities_Ha[34]};

// From the GUI
int_lines[] = {int_lines[],183,218,197,241,249,205,227};
int_points[] = {int_points[],83,95,77,89};

ext_lines[] = {ext_lines[],201,293,275,271,231,178,187,179,195,285,173,263,196,284};
ext_points[] = {ext_points[],79,107,101,73};

outlet_lines[0] = 240;

//-----------------------------------------------------------------------------
//Thrid extrusion along the direction of the B field (wall);

extruded_entities_Ha[] = Extrude {0.0,t,0.0} {Surface{extruded_entities_Ha[0],extruded_entities_Ha[6],extruded_entities_Ha[12],extruded_entities_Ha[18],extruded_entities_Ha[24],extruded_entities_Ha[30],15,27,39};
									Layers{N_w}; Recombine;};

//For i In {0:53}
//	Printf("%f",extruded_entities_Ha[i]);
//EndFor

//Physical tags from the extruded list

wall_volumes[] = {wall_volumes[],
		extruded_entities_Ha[1],extruded_entities_Ha[7],extruded_entities_Ha[13],extruded_entities_Ha[25],extruded_entities_Ha[31],extruded_entities_Ha[37],extruded_entities_Ha[43],extruded_entities_Ha[49]}; 
fluid_volumes[1] = extruded_entities_Ha[19];

int_surfaces[] = {int_surfaces[],
		extruded_entities_Ha[12],extruded_entities_Ha[20], extruded_entities_Ha[22], extruded_entities_Ha[23], extruded_entities_Ha[42]};
outlet_surfaces[1] = extruded_entities_Ha[21];
ext_surfaces[] = {ext_surfaces[],
		extruded_entities_Ha[2],extruded_entities_Ha[8], extruded_entities_Ha[9],extruded_entities_Ha[28], extruded_entities_Ha[33],extruded_entities_Ha[34], extruded_entities_Ha[38], extruded_entities_Ha[41], extruded_entities_Ha[47], extruded_entities_Ha[52], extruded_entities_Ha[53]};

// From the GUI
int_lines[] = {int_lines[],337,381,329,373,439,462,461,307,351,315,359,350,352};
int_points[] = {int_points[],123,139,129,133,113,117,165,175};

ext_lines[] = {ext_lines[],333,425,327,417,403,395,407,311,305,310,437,483,451,440,473,484,442,495,328,416};
ext_points[] = {ext_points[],155,119,108,109,149,145,156,185};

outlet_lines[1] = 372;

//--------------------------------------------------------------------------
//Fourth extrusion along the direction of the B field (inlet channel volume);

extruded_entities_Ha[] = Extrude {0.0,2/r,0.0} {Surface{extruded_entities_Ha[0],extruded_entities_Ha[6],extruded_entities_Ha[12],extruded_entities_Ha[18],extruded_entities_Ha[24],extruded_entities_Ha[30],extruded_entities_Ha[36],extruded_entities_Ha[42],extruded_entities_Ha[48]};
										Layers{one[],layer[]}; Recombine;};
										
//Physical tags from the extruded list
wall_volumes[] = {wall_volumes[],
		extruded_entities_Ha[1],extruded_entities_Ha[7],extruded_entities_Ha[25],extruded_entities_Ha[31],extruded_entities_Ha[37],extruded_entities_Ha[49]}; 
fluid_volumes[] = {fluid_volumes[], extruded_entities_Ha[13],extruded_entities_Ha[19],extruded_entities_Ha[43]};

int_surfaces[] = {int_surfaces[],
		extruded_entities_Ha[4],extruded_entities_Ha[10],extruded_entities_Ha[12], extruded_entities_Ha[26], extruded_entities_Ha[32], extruded_entities_Ha[42], extruded_entities_Ha[44], extruded_entities_Ha[46]};
outlet_surfaces[2] = extruded_entities_Ha[21];
inlet_surfaces[0] = extruded_entities_Ha[47];
ext_surfaces[] = {ext_surfaces[],
		extruded_entities_Ha[2],extruded_entities_Ha[8], extruded_entities_Ha[9],extruded_entities_Ha[28], extruded_entities_Ha[33],extruded_entities_Ha[34], extruded_entities_Ha[38], extruded_entities_Ha[41], extruded_entities_Ha[52], extruded_entities_Ha[53]};

// From the GUI
int_lines[] = {int_lines[],535,579,527,571,550,549,505,649,660,671,637,659,517,513,561,557,548};
int_points[] = {int_points[],201,217,207,211,191,195,243,253};

ext_lines[] = {ext_lines[],526,614,525,615,503,593,635,681,682,638,640,693,623,531,605,601,509,508};
ext_points[] = {ext_points[],197,233,223,227,187,186,234,263};

outlet_lines[2] = 570;

//-----------------------------------------------------------------------------
//Fith extrusion along the direction of the B field (wall);
extruded_entities_Ha[] = Extrude {0.0,t,0.0} {Surface{extruded_entities_Ha[0],extruded_entities_Ha[6],extruded_entities_Ha[12],extruded_entities_Ha[18],extruded_entities_Ha[24],extruded_entities_Ha[30],extruded_entities_Ha[36],extruded_entities_Ha[42],extruded_entities_Ha[48]};
									Layers{N_w}; Recombine;};

//Physical tags from the extruded list
wall_volumes[] = {wall_volumes[],
		extruded_entities_Ha[1],extruded_entities_Ha[7], extruded_entities_Ha[13],extruded_entities_Ha[25],extruded_entities_Ha[31],extruded_entities_Ha[37],extruded_entities_Ha[43],extruded_entities_Ha[49]}; 
fluid_volumes[] = {fluid_volumes[], extruded_entities_Ha[19]};

int_surfaces[] = {int_surfaces[],
		extruded_entities_Ha[10], extruded_entities_Ha[15], extruded_entities_Ha[32]};
outlet_surfaces[3] = extruded_entities_Ha[21];
ext_surfaces[] = {ext_surfaces[],
		extruded_entities_Ha[2],extruded_entities_Ha[8], extruded_entities_Ha[9],extruded_entities_Ha[28], extruded_entities_Ha[33],extruded_entities_Ha[34], extruded_entities_Ha[36], extruded_entities_Ha[38], extruded_entities_Ha[41], extruded_entities_Ha[42], extruded_entities_Ha[47], extruded_entities_Ha[48], extruded_entities_Ha[52], extruded_entities_Ha[53]};

// From the GUI
int_lines[] = {int_lines[],733,777,725,769,711,755,746};
int_points[] = {int_points[],279,295,269,285};

ext_lines[] = {ext_lines[],723,813,724,812,791,748,701,833,835,857,879,858,836,880,838,847,869,891,835,748,729,803,799,821,704,792};
ext_points[] = {ext_points[],275,311,723,264,265,273,289,301,305,312,321,331,341};

outlet_lines[3] = 768;


//--------------------------------------------------------
//Sixth extrusion along the direction of the B field (PbLi volume);
extruded_entities_Ha[] = Extrude {0.0,(r-1)/r-t,0.0} {Surface{extruded_entities_Ha[0],extruded_entities_Ha[6],extruded_entities_Ha[12],extruded_entities_Ha[18],extruded_entities_Ha[24],extruded_entities_Ha[30]};
									Layers{one[],layer[]}; Recombine;};
									
//Physical tags from the extruded list
wall_volumes[] = {wall_volumes[],
		extruded_entities_Ha[1],extruded_entities_Ha[7], extruded_entities_Ha[13],extruded_entities_Ha[25],extruded_entities_Ha[31]}; 
fluid_volumes[] = {fluid_volumes[], extruded_entities_Ha[19]};

int_surfaces[] = {int_surfaces[],
		extruded_entities_Ha[10], extruded_entities_Ha[15], extruded_entities_Ha[18], extruded_entities_Ha[32]};
outlet_surfaces[4] = extruded_entities_Ha[21];
ext_surfaces[] = {ext_surfaces[],
		extruded_entities_Ha[2],extruded_entities_Ha[5],extruded_entities_Ha[8], extruded_entities_Ha[9],extruded_entities_Ha[17],extruded_entities_Ha[28], extruded_entities_Ha[29], extruded_entities_Ha[33],extruded_entities_Ha[34]};
							
// From the GUI
int_lines[] = {int_lines[],931,975,923,967,909,953,944,966};
int_points[] = {int_points[],357,373,347,363};	

ext_lines[] = {ext_lines[],922,921,1010,1011,946,899,902,905,904,913,957,1001,997,927,1019,990,989};
ext_points[] = {ext_points[],353,343,342,351,367,383,379};	

			
									
//-----------------------------------------------------------------------------
//Seventh extrusion along the direction of the B field (PbLi volume);
extruded_entities_wall[] = Extrude {0.0,t,0.0} {Surface{extruded_entities_Ha[0],extruded_entities_Ha[6],extruded_entities_Ha[12],extruded_entities_Ha[18],extruded_entities_Ha[24],extruded_entities_Ha[30]};
									Layers{N_w}; Recombine;};

//Physical tags from the extruded list
wall_volumes[] = {wall_volumes[],
		extruded_entities_wall[1],extruded_entities_wall[7], extruded_entities_wall[13], extruded_entities_wall[19], extruded_entities_wall[25],extruded_entities_wall[31]}; 

ext_surfaces[] = {ext_surfaces[],
		extruded_entities_wall[0],extruded_entities_wall[2],extruded_entities_wall[5],extruded_entities_wall[6],extruded_entities_wall[8],extruded_entities_wall[9], extruded_entities_wall[12],extruded_entities_wall[17], extruded_entities_wall[18], extruded_entities_wall[21], extruded_entities_wall[24], extruded_entities_wall[28],extruded_entities_wall[29], extruded_entities_wall[30],extruded_entities_wall[33],extruded_entities_wall[34]};
				
ext_lines[] = {ext_lines[],1098,1063,1054,1059,1053,1055,1031,1037,1036,1045,1034,1033,1078,1076,1077,1089,1133,1129,1121,1122,1120,1099,1143,1011,1142,1151,1107,1032};
ext_points[] = {ext_points[],401,405,421,437,390,391,399,395,415,431,411,427,421,437};				
							
//-------------Definition of physical tags from the lists---------------

Physical Point ("exterior") = {ext_points[]};//
Physical Line ("exterior") = {ext_lines[]};
Physical Surface ("exterior") = {ext_surfaces[]};

Physical Point ("interior") = {int_points[]};
Physical Line ("interior") = {int_lines[]};
Physical Surface ("interior") = {int_surfaces[]};

Physical Line ("outlet") = {outlet_lines[]};
Physical Surface ("outlet") = {outlet_surfaces[]};

Physical Surface ("inlet") = {inlet_surfaces[]};

Physical Volume("fluid") = {fluid_volumes[]};
Physical Volume("wall") = {wall_volumes[]};
