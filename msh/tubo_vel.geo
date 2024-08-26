//
// Tube with normal velocity 
//

// Element size
lc = 0.01;

// Corners
Point(1) = { 0, 0,     0, lc};
Point(2) = { 1, 0,     0, lc};
Point(3) = { 1, 0.02,  0, lc};
Point(4) = { 0, 0.02,  0, lc};

// Edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,Ar,1,1.028,340.0,400.0") = {1};

// Normal velocity (left)
// closed-open
// f1 = 85 Hz
//Physical Curve("Vn,0.0,160.0,0.0") = {4};

// Atenuation
Physical Curve("Cn,1.0") = {2,4};

// Boundary conditions - Open
// Physical Curve("Open") = {2};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "tubo_vel.msh";
