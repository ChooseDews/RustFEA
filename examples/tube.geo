// Gmsh project created on Sat Sep 28 00:03:43 2024
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 1, 0, 2*Pi};
//+
Transfinite Curve {2, 1} = 20 Using Progression 1;

//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
//+
Extrude {0, 0, 20} {
  Surface{1}; Layers {20}; Recombine;
}
//+
Physical Surface("end1", 7) = {4};
//+
Physical Surface("end2", 8) = {1};
//+
Physical Curve("end_outer", 9) = {4};
//+
Physical Curve("end_inner", 10) = {6};
