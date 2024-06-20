
// (Externally written) Mesh parameters^M
r=1.2626048580501072;
R=1.0;
N_r=12;
n=4;
N_a=40;
N_L=12;
L=4.0;
p=0.9;
q=0.5;


// Circle

//+Centre
Point(1) = {0.0, 0.0, -0.5*L, 1.0};

//+
Point(2) = {0.0, R, -0.5*L, 1.0};
//+
Point(3) = {R, 0.0, -0.5*L, 1.0};
//+
Point(4) = {0.0, -R, -0.5*L, 1.0};
//+
Point(5) = {-R, 0.0, -0.5*L, 1.0};

//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};

// internal square

//+
Point(6) = {0.0, q, -0.5*L, 1.0};
//+
Point(7) = {q, 0.0, -0.5*L, 1.0};
//+
Point(8) = {0.0, -q, -0.5*L, 1.0};
//+
Point(9) = {-q, 0.0, -0.5*L, 1.0};

//Intermediate circle

//+
Point(10) = {0.0, p, -0.5*L, 1.0};
//+
Point(11) = {p, 0.0, -0.5*L, 1.0};
//+
Point(12) = {0.0, -p, -0.5*L, 1.0};
//+
Point(13) = {-p, 0.0, -0.5*L, 1.0};

// sides of the square

//+
Line (5) = {6,7};
//+
Line (6) = {7,8};
//+
Line (7) = {8,9};
//+
Line (8) = {9,6};

//Radious between internal circle and square

//+
Line (9) = {10,6};
//+
Line (10) = {11,7};
//+
Line (11) = {12,8};
//+
Line (12) = {13,9};

// BL circle

Circle(13) = {10, 1, 11};
//+
Circle(14) = {11, 1, 12};
//+
Circle(15) = {12, 1, 13};
//+
Circle(16) = {13, 1, 10};

//Radious in the BL

//+
Line (17) = {2,10};
//+
Line (18) = {3,11};
//+
Line (19) = {4,12};
//+
Line (20) = {5,13};

//Surfaces

//+
Curve Loop(1) = {9, 5, -10, -13}; 		Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 6, -11, -14}; 		Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 7, -12, -15}; 		Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 8, -9, -16}; 		Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 6, 7, 8}; 			Plane Surface(5) = {5};
//+
Curve Loop(6) = {17, 13, -18, -1}; 		Plane Surface(6) = {6};
//+
Curve Loop(7) = {18, 14, -19, -2}; 		Plane Surface(7) = {7};
//+
Curve Loop(8) = {19, 15, -20, -3}; 		Plane Surface(8) = {8};
//+
Curve Loop(9) = {20, 16, -17, -4}; 		Plane Surface(9) = {9};
//Azimutal nodes

Transfinite Curve {1, 2, 3, 4, 5, 6 , 7, 8, 13, 14, 15, 16} = N_a/4;

//Radial nodes (internal)

n = (p-q)/(R-p)*(r^(-N_r)-1)/(1-r);

Transfinite Curve {9, 10, 11, 12} = n;

//Radial nodes (BL)

Transfinite Curve {17, 18, 19, 20} = N_r Using Progression r;

// Superficial Mesh

Transfinite Surface {:};
//Recombine Surface {1, 2, 3, 4, 5};

// Extrusion
Extrude {0.0,0.0,L} {Surface{:}; Layers{N_L};} 

//Physical tags
//+
Physical Surface("inlet") = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Physical Line("inlet") = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 ,18 ,19, 20};
Physical Point("inlet") = {1, 6, 7, 8, 9, 10, 11, 12, 13};
//+
Physical Surface("outlet") = {152, 218, 42, 86, 64, 174, 196, 108, 130};
Physical Line("outlet") = {23, 45, 67, 89,22, 24, 46, 68, 47, 69, 91, 25, 132, 134, 156, 178};
Physical Point("outlet") = {27, 23, 37, 47, 14, 15, 19, 33, 43};
//+
Physical Surface("wall") = {195, 173, 217, 151};
Physical Line("wall") = {1, 2, 3, 4, 157, 179, 201, 135, 137, 168, 190, 146};
Physical Point("wall") = {2, 3, 4, 5, 58, 69, 80, 48};
//+
Physical Volume("Fluid") = {1, 2, 3, 4, 5, 6, 7, 8, 9};

