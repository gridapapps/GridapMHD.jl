//Geometry parameters

//Radious
R = 1.0;


//Pipe lenght
L = 2.5; 

//side of the internal square
p = 0.5;


// Circle

//+Centre
Point(1) = {0.0, 0.0, -L, 1.0};

//+
Point(2) = {0.0, R, -L, 1.0};
//+
Point(3) = {R, 0.0, -L, 1.0};
//+
Point(4) = {0.0, -R, -L, 1.0};
//+
Point(5) = {-R, 0.0, -L, 1.0};

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
Point(6) = {0.0, p, -L, 1.0};
//+
Point(7) = {p, 0.0, -L, 1.0};
//+
Point(8) = {0.0, -p, -L, 1.0};
//+
Point(9) = {-p, 0.0, -L, 1.0};

// sides of the square

//+
Line (5) = {6,7};
//+
Line (6) = {7,8};
//+
Line (7) = {8,9};
//+
Line (8) = {9,6};

//Radious between circle and square

//+
Line (9) = {2,6};
//+
Line (10) = {3,7};
//+
Line (11) = {4,8};
//+
Line (12) = {5,9};


//Surfaces

//+
Curve Loop(1) = {9, 5, -10, -1}; 		Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 6, -11, -2}; 		Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 7, -12, -3}; 		Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 8, -9, -4}; 		Plane Surface(4) = {4};
//+
Curve Loop(5) = {5, 6, 7, 8}; 			Plane Surface(5) = {5};

//Azimutal nodes

Transfinite Curve {1, 2, 3, 4, 5, 6 , 7, 8} = 10;

//Radial nodes

Transfinite Curve {9, 10, 11, 12} = 25 Using Progression 1.25;

// Superficial Mesh

Transfinite Surface {1, 2, 3, 4, 5};
//Recombine Surface {1, 2, 3, 4, 5};

// Extrusion
Extrude {0.0,0.0,L} {Surface{1, 2, 3, 4, 5}; Layers{10};} 

//Physical tags
//+
Physical Surface("inlet") = {1, 2, 3, 4, 5};
Physical Line("inlet") = {5, 6, 7, 8, 9, 10, 11, 12};
Physical Point("inlet") = {6, 7, 8, 9};
//+
Physical Surface("outlet") = {100, 34, 56, 78, 122};
Physical Line("outlet") = {81, 15, 37, 59, 14, 16, 38, 60};
Physical Point("outlet") = {39, 29, 15, 11};
//+
Physical Surface("wall") = {99, 33, 55, 77};
Physical Line("wall") = {1, 2, 3, 4, 17, 39, 61, 83, 19, 28, 50, 72};
Physical Point("wall") = {2, 43, 33, 19, 5, 4, 3, 10};
//+
Physical Volume("Fluid") = {1, 2, 3, 4, 5};

