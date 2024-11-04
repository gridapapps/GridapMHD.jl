//Geometry parameters

// (Externally written) Mesh parameters
r=1.2626048580501072;
R=1.0;
N_r=12;
n=4;
n_c=0.2;
N_a=40;
N_L=12;
L=4.0;
p=0.9;
q=0.5;


// Circle

//+Centre
Point(1) = {0.0, 0.0, -0.5*L, 1.0};

//4 poits per shell

//+
Point(2) = {0.0, R, -0.5*L, 1.0};
//+
Point(3) = {R, 0.0, -0.5*L, 1.0};
//+
Point(4) = {0.0, -R, -0.5*L, 1.0};
//+
Point(5) = {-R, 0.0, -0.5*L, 1.0};

//+
Point(6) = {0.0, p*R, -0.5*L, 1.0};
//+
Point(7) = {p*R, 0.0, -0.5*L, 1.0};
//+
Point(8) = {0.0, -p*R, -0.5*L, 1.0};
//+
Point(9) = {-p*R, 0.0, -0.5*L, 1.0};

//+
Point(10) = {0.0, q*R, -0.5*L, n_c*R};
//+
Point(11) = {q*R, 0.0, -0.5*L, n_c*R};
//+
Point(12) = {0.0, -q*R, -0.5*L, n_c*R};
//+
Point(13) = {-q*R, 0.0, -0.5*L, n_c*R};

//Curves

//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};

//+
Circle(5) = {6, 1, 7};
//+
Circle(6) = {7, 1, 8};
//+
Circle(7) = {8, 1, 9};
//+
Circle(8) = {9, 1, 6};

//+
Circle(9) = {10, 1, 11};
//+
Circle(10) = {11, 1, 12};
//+
Circle(11) = {12, 1, 13};
//+
Circle(12) = {13, 1, 10};

//Radious for the transfinite cells

Line(13) = {2,6};
Line(14) = {3,7};
Line(15) = {4,8};
Line(16) = {5,9};

//Extra shell for smother transition

t = 0.5*(p+q);
n_t = 2*3.141592*R/N_a;

//+
Point(14) = {0.0, t*R, -0.5*L, n_t*R};
//+
Point(15) = {t*R, 0.0, -0.5*L, n_t*R};
//+
Point(16) = {0.0, -t*R, -0.5*L, n_t*R};
//+
Point(17) = {-t*R, 0.0, -0.5*L, n_t*R};

//+
Circle(18) = {14, 1, 15};
//+
Circle(19) = {15, 1, 16};
//+
Circle(20) = {16, 1, 17};
//+
Circle(21) = {17, 1, 14};

//Surfaces

//+
Curve Loop(1) = {1, 14,-5,-13};			 		Plane Surface(1) = {1};
Curve Loop(2) = {2, 15,-6,-14};			 		Plane Surface(2) = {2};
Curve Loop(3) = {3, 16,-7,-15};			 		Plane Surface(3) = {3};
Curve Loop(4) = {4, 13,-8,-16};			 		Plane Surface(4) = {4};
Curve Loop(5) = {5, 6, 7, 8, 18, 19, 20, 21};	Plane Surface(5) = {5};
Curve Loop(6) = {18, 19, 20, 21, 9, 10, 11, 12};Plane Surface(6) = {6};
Curve Loop(7) = {9, 10, 11, 12};				Plane Surface(7) = {7};

//Recombine Surface{:};

Transfinite Curve {13,14,15,16} = N_r Using Progression r;
Transfinite Curve {1,2,3,4,5,6,7,8} = N_a/4;

Transfinite Surface {1,2,3,4}; 

//+ Extrusion
//Extrude {0.0,0.0,L}{ Surface{:}; Layers{N_L}; Recombine; } 
Extrude {0.0,0.0,L}{ Surface{:}; Layers{N_L}; }

//Physical tags
//+
Physical Surface("inlet") = {1, 2, 3, 4, 5, 6, 7};
Physical Line("inlet") = {5, 6, 7, 8, 9, 10, 11, 12, 18, 19, 20, 21, 13, 14, 15, 16};
Physical Point("inlet") = {1, 10, 14, 6, 11, 15, 7, 12, 16, 8, 13, 17, 9};
//+
Physical Surface("outlet") = {215, 87, 151, 43, 109, 193, 65};
Physical Line("outlet") = {25, 47, 69, 91, 115, 116, 69, 118, 157, 158, 159, 160, 26, 24, 46, 68, 117};
Physical Point("outlet") = {19, 97, 64, 29, 99, 66, 24, 104, 71, 36, 109, 76, 43};
//+
Physical Surface("wall") = {74, 52, 30, 96};
Physical Line("wall") = {1, 2, 3, 4, 23, 45, 67, 89, 28, 29, 51, 73};
Physical Point("wall") = {2, 3, 4, 5, 18, 20, 32, 43};
//+
Physical Volume("Fluid") = {1,2,3,4,5,6,7};
