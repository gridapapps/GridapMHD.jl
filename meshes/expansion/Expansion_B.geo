// Gmsh project created on Tue Jan 16 14:55:13 2023

// (Externally written) Mesh parameters
R_Ha=1.2972667023986404;
R_side=1.2432187557741343;
N_Ha=20;
n_Ha=4;
N_side=20;
n_side=5;
n_inlet=20;
n_outlet=40;
r_exp=1.075;
r=4.0;
betta=0.25;
L_in=1.0;
L_out=2.0;




//2D points in z=-betta and then strude

//+
Point(1) = {-L_in, 0, -betta, 1.0};
//+
Point(2) = {-L_in, 1/r, -betta, 1.0};
//+
Point(3) = {0, 1/r, -betta, 1.0};
//+
Point(4) = {0, (1+1/r)*0.5, -betta, 1.0};
//+
Point(5) = {0, 1.0, -betta, 1.0};
//+
Point(6) = {L_out, 1.0, -betta, 1.0};
//+
Point(7) = {L_out, (1+1/r)*0.5, -betta, 1.0};
//+
Point(8) = {L_out, 1/r, -betta, 1.0};
//+
Point(9) = {L_out, 0, -betta, 1.0};
//+
Point(10) = {L_out, -1/r, -betta, 1.0};
//+
Point(11) = {L_out, -(1+1/r)*0.5, -betta, 1.0};
//+
Point(12) = {L_out, -1.0, -betta, 1.0};
//+
Point(13) = {0, -1.0, -betta, 1.0};
//+
Point(14) = {0, -(1+1/r)*0.5, -betta, 1.0};
//+
Point(15) = {0, -1/r, -betta, 1.0};
//+
Point(16) = {-L_in, -1/r, -betta, 1.0};

//Extra point in the center of the domain
//+
Point(17) = {0, 0, -betta, 1.0};

//Lines joining points

//+
Line(1) = {1, 2};
//+
Line(2) = {2,3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 11};
//+
Line(11) = {11, 12};
//+
Line(12) = {12, 13};
//+
Line(13) = {13, 14};
//+
Line(14) = {14, 15};
//+
Line(15) = {15, 16};
//+
Line(16) = {16, 1};
//+
Line(17) = {3, 17};
//+
Line(18) = {17, 15};


//2D mesh in the z = -betta plane

//+
Transfinite Curve {-2, 15} = n_inlet Using Progression r_exp; //Axial cells in the inlet channel
//+
Transfinite Curve {5, -12} = n_outlet Using Progression r_exp; //Axial cells in the outlet channel
//+
Transfinite Curve {-1, 17, 8, 16, -18, -9} = N_Ha/2 Using Progression R_Ha; //Cells along the B direction in the inlet channel (x2)
//+
Transfinite Curve {3, -7, 13, -11, -14, 10, -4, 6} = N_Ha/2 Using Progression R_Ha; //Cells along the B direction in the outle channel (x4)

//Surface definition

//+
Curve Loop(1) = {1, 2, 17, 18, 15, 16};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {-17, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, -18};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {1} = {2, 3, 15, 16};
//+
Transfinite Surface {2} = {5, 6, 12, 13};
//+
Recombine Surface {1,2}; //This recombines the triangles into rectangles
//+

//Extrude using also a geometric progression along the direction perpendicular to B

N = N_side/2;
one[0] = 1;
layer[0] = 0.5*(1 - R_side) / (1-R_side^N);

For i In {1:N-1}
   one[i] = 1;
   layer[i] = layer[i-1] + layer[0] * R_side^i;
EndFor

For i In {0:N-1}
   one[i+N] = 1;
   layer[i+N] = layer[i+N-1] + layer[0]*R_side^(N-i-1);
EndFor

Extrude {0.0,0.0,2*betta} {
    Surface{1,2}; Layers{one[],layer[]}; Recombine; 
}


//Tag definition for gridap

//+
Physical Surface("inlet") = {29, 49};
//+
Physical Surface("outlet") = {85, 89, 93, 97, 101, 105};
//+
Physical Surface("wall") = {1, 2, 33, 45, 109, 81, 50, 122, 77, 73, 117, 113};
//+
Physical Volume("PbLi") = {1, 2};
//+
Physical Line("inlet") = {27};
Physical Line("outlet") = {84,88,92,96,100};
Physical Line("wall") = {1,16,20,25,28,44,2,15,21,24,3,4,17,18,13,14,53,54,22,23,63,64,76,72,32,40,112,108,5,12,55,62,6,7,8,9,10,11,56,57,58,59,60,61,80,104};
//+
Physical Point("wall") = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,27,31,35,41,45,49,53,57,61,65,69,73,77,81};
