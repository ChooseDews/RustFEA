// Gmsh project created on Thu Oct 10 23:06:15 2024
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {-2, -1.1117236434354394, -0.025, 1, 2, 0};
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
//+
Transfinite Surface {2} = {13, 11, 10, 12};
//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 2, 3, 4} = 1 Using Progression 1;
//+
Extrude {0, 0, 0.25} {
  Surface{2}; Layers {1}; Recombine;
}
//+
Extrude {0, 0, 0.3} {
  Surface{1}; Layers {1}; Recombine;
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
//+
//+
Transfinite Curve {8} = 30 Using Bump 15;
//+
Transfinite Curve {5} = 30 Using Bump 1;
//+
Physical Surface("back", 31) = {7};
//+
Transfinite Curve {6, 7} = 12 Using Progression 0.7;
//+
Physical Curve("back_line", 32) = {16};
//+
Transfinite Surface {9} = {20, 19, 2, 3};
//+
Transfinite Volume{2} = {19, 20, 21, 18, 1, 4, 3, 2};
//+
MeshSize {20, 3, 2, 19} = 0.05;
//+
Transfinite Curve {21, 20, 2, 18} = 10 Using Progression 1;
//+
Transfinite Curve {4, 22, 24, 17} = 10 Using Progression 1;
//+
Transfinite Curve {21, 20, 2, 18, 19, 1, 23, 3, 24, 4, 17, 22} = 20 Using Progression 1;
//+
Transfinite Curve {20, 22, 18, 17} = 100 Using Progression 1;
//+
Transfinite Curve {2, 20, 23, 3, 22, 24, 4, 17, 19, 1, 18, 21} = 10 Using Progression 1;
//+
Recombine Surface {12, 1, 8};
//+
Transfinite Curve {20, 22, 18, 2, 21, 23, 3, 4, 24, 19, 17, 1} = 20 Using Progression 1;
