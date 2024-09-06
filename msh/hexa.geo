Merge "geom3.brep";
Transfinite Line "*" = 20 Using Bump 1.;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

// Material
Physical Volume("Material,Ar,1,1.028,340.0,400.0") = {1,2,3};

// Normal velocity (left)
Physical Surface("Vn,1E-3,340.0,0.0") = {1};

// Probe
Physical Curve("Probe") = {9};

// Atenuation
Physical Surface("Yn,1.0") = {1,14};


