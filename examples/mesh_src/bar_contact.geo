// Gmsh project created on Tue Oct  1 23:03:05 2024
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 0.5, 0};
//+
Transfinite Surface {1};
//+
Extrude {0, 0, 10} {
  Surface{1}; Layers {27}; Recombine;
}
//+
Physical Surface("free_end", 13) = {6};
//+
Physical Surface("fixed_end", 14) = {1};
//+
//+
Rectangle(7) = {-1.01, -0.25, 2, 1, 1, 0};
//+
Transfinite Surface {7} = {11, 10, 9, 12};
//+
Extrude {0, 0, 1} {
  Surface{7}; Layers {1}; Recombine;
}
//+
Transfinite Curve {13, 16, 15, 14} = 1 Using Progression 1;
//+
Physical Volume("bar_body", 25) = {1};
//+
Physical Volume("block_body", 26) = {2};
//+
Physical Surface("bar_bottom", 27) = {5};
//+
Physical Surface("block_top", 28) = {9};
