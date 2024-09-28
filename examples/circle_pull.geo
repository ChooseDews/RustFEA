// Gmsh project created on Sat Sep 28 01:16:18 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 1, 0, Pi/2};
//+
Circle(2) = {0, 0, 0, 0.25, 0, Pi/2};
//+
Line(3) = {1, 3};
//+
Line(4) = {4, 2};
//+
Curve Loop(1) = {1, -4, -2, -3};
//+
Plane Surface(1) = {1};//+
Transfinite Surface {1};
//+
Transfinite Curve {1, 2} = 30 Using Progression 1;
//+
Transfinite Curve {4, 3} = 25 Using Progression 1;
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers {10}; Recombine;
}

//+
Physical Surface("side1", 13) = {5};
//+
Physical Surface("side2", 14) = {3};
//+
Physical Surface("inside", 15) = {4};
//+
Physical Surface("outside", 16) = {2};
