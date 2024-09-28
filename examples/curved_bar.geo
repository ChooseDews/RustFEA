// Gmsh project created on Fri Sep 27 17:42:37 2024
SetFactory("OpenCASCADE");
//+
Torus(1) = {0, 0, 0, 0.5, 0.2, 2*Pi};
//+
Box(2) = {0, -2, -0.5, 1, 4, 1};
//+
BooleanIntersection{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+remove another half to leave 1/4 of the torus
Box(3) = {0, 0, -0.5, 1, 4, 1};
//+
BooleanIntersection{ Volume{1}; Delete; }{ Volume{3}; Delete; }
//+
Physical Surface("end1", 4) = {3};
//+
Physical Surface("end2", 5) = {2};
//+
Physical Surface("surface", 6) = {1};
