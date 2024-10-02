// Gmsh project created on Tue Oct  1 23:03:05 2024
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 0.5, 0};
//+
Transfinite Surface {1};
//+
Extrude {0, 0, 10} {
  Surface{1}; Layers {10}; Recombine;
}
//+
Physical Surface("free_end", 13) = {6};
//+
Physical Surface("fixed_end", 14) = {1};
