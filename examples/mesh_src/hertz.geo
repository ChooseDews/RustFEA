// Gmsh project created on Thu Oct 10 23:06:15 2024
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-2, -1, 0, 1, 2, 0};
Point(9) = {0, 0, 0, 1.0};
//+
Point(10) = {-0.2, -1, 0, 1.0};
Point(11) = {-0.2, 1, 0, 1.0};
//+
Line(5) = {10, 11};
//+
//+
Point(12) = {-0.70710678118, -0.70710678118, 0, 1.0};
//+
Point(13) = {-0.70710678118, 0.70710678118, 0, 1.0};
//+
Line(6) = {10, 12};
//+
Line(7) = {11, 13};
//+
Circle(8) = {12, 9, 13};
//+
Curve Loop(2) = {6, 8, -7, -5};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {8, 5} = 30 Using Progression 1;
//+
Transfinite Curve {6, 7} = 30 Using Progression 1;
//+
Transfinite Surface {2} = {13, 11, 10, 12};
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 2, 3, 4} = 1 Using Progression 1;
//+
Extrude {0, 0, 0.25} {
  Surface{2}; Surface{1}; Layers {1}; Recombine;
}
//+
Physical Surface("top", 25) = {6};
//+
Physical Surface("bottom", 26) = {11};
//+
Physical Surface("contact_a", 27) = {4};
//+
Physical Surface("contact_b", 28) = {9};
//+
Physical Volume("block", 29) = {2};
//+
Physical Volume("cylinder", 30) = {1};
