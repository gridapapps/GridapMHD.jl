//Geometry parameters

//Radious
R = 1.0;		//Radious
p =0.95;		//BL thickness
q =0.7;		//Core region

L = 4.0;		//Pipe lenght 

n_c = 0.1;	//Core cell size
n_BL = 0.01;	//BL cell size
n_L = 12;	//Axial cells


// Circle

//+Centre
Point(1) = {0.0, 0.0, -0.5*L, 1.0};

//4 poits per shell

//+
Point(2) = {0.0, R, -0.5*L, n_BL*R};
//+
Point(3) = {R, 0.0, -0.5*L, n_BL*R};
//+
Point(4) = {0.0, -R, -0.5*L, n_BL*R};
//+
Point(5) = {-R, 0.0, -0.5*L, n_BL*R};

//+
Point(6) = {0.0, p*R, -0.5*L, n_BL*R};
//+
Point(7) = {p*R, 0.0, -0.5*L, n_BL*R};
//+
Point(8) = {0.0, -p*R, -0.5*L, n_BL*R};
//+
Point(9) = {-p*R, 0.0, -0.5*L, n_BL*R};

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


//Surfaces

//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8}; 		Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8, 9, 10, 11, 12};	Plane Surface(2) = {2};
Curve Loop(3) = {9, 10, 11, 12};				Plane Surface(3) = {3};

//Recombine Surface{1,2,3};

//+ Extrusion
//Extrude {0.0,0.0,L}{ Surface{:}; Layers{n_L}; Recombine; }
Extrude {0.0,0.0,L}{ Surface{:}; Layers{n_L}; }

//Physical tags
//+
Physical Surface("inlet") = {1, 2, 3};
Physical Line("inlet") = {5, 6, 7, 8, 9, 10, 11, 12};
Physical Point("inlet") = {1, 6, 7, 8, 9, 10, 11, 12, 13};
//+
Physical Surface("outlet") = {54, 96, 118};
Physical Line("outlet") = {18, 19, 20, 21, 60, 61, 62, 63};
Physical Point("outlet") = {15, 34, 36, 41, 46, 67, 69, 74, 79};
//+
Physical Surface("wall") = {25, 29, 33, 37};
Physical Line("wall") = {1, 2, 3, 4, 14, 15, 16, 17, 23, 24, 28, 32};
Physical Point("wall") = {2, 3, 4, 5, 14, 16, 21, 26};
//+
Physical Volume("Fluid") = {1,2,3};

