// Gmsh project created on Wed Jan 19 14:55:13 2022

//Puntos de la figura 2D  en z=-1.0 para luego extruir

//+
Point(1) = {-8, 0, -1.0, 1.0};
//+
Point(2) = {-8, 0.25, -1.0, 1.0};
//+
Point(3) = {0, 0.25, -1.0, 1.0};
//+
Point(4) = {0, 0.625, -1.0, 1.0};
//+
Point(5) = {0, 1.0, -1.0, 1.0};
//+
Point(6) = {8, 1.0, -1.0, 1.0};
//+
Point(7) = {8, 0.625, -1.0, 1.0};
//+
Point(8) = {8, 0.25, -1.0, 1.0};
//+
Point(9) = {8, 0, -1.0, 1.0};
//+
Point(10) = {8, -0.25, -1.0, 1.0};
//+
Point(11) = {8, -0.625, -1.0, 1.0};
//+
Point(12) = {8, -1.0, -1.0, 1.0};
//+
Point(13) = {0, -1.0, -1.0, 1.0};
//+
Point(14) = {0, -0.625, -1.0, 1.0};
//+
Point(15) = {0, -0.25, -1.0, 1.0};
//+
Point(16) = {-8, -0.25, -1.0, 1.0};

//Punto extra en el centro
//+
Point(17) = {0, 0, -1.0, 1.0};

//Lineas que unen los puntos

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

//Dos líneas extra para poder definir luego las superficies con sus hexahedros
//+
Line(17) = {3, 17};
//+
Line(18) = {17, 15};


//Defino los nodos de las líneas

//+
Transfinite Curve {-2, 15, 5, -12} = 100 Using Progression 1.02; //Nodos axiales
//+
Transfinite Curve {-1, 17, 8, 16, -18, -9} = 30 Using Progression 1.4; //Nodos dirección Hartmann en el canal pequeño
//+
Transfinite Curve {3, -7, 13, -11, -14, 10, -4, 6} = 30 Using Progression 1.4; //Nodos dirección Hartmann en el canal pequeño

//Defino la superfice

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
Recombine Surface {1,2}; //Esto hace desaparecer los triangulos
//+
Extrude {0.0,0.0,2.0} {
    Surface{1,2}; Layers{{8,23,8}, {0.03,0.97,1}}; Recombine; //Pongo 8 nodos en las capas limites side para Ha = 1000
}

//Defino las tags para las condiciones de contorno (los números de las superficies los genera gmsh automáticamente al extruir)

//+
Physical Surface("inlet") = {29, 49};
//+
Physical Surface("outlet") = {85, 89, 93, 97, 101, 105};
//+
Physical Surface("wall") = {1, 2, 33, 45, 109, 81, 50, 122, 77, 73, 117, 113};
//+
Physical Volume("PbLi") = {1, 2};
