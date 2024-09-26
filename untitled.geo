//+
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 1, 0, 0, 0.02, 2*Pi};
//+
Cylinder(2) = {0, 0, 0, 1, 0, 0, 0.02, 2*Pi};
//+
Physical Volume("Material,ar,1,1.028,340.0,400.0") = {1};
