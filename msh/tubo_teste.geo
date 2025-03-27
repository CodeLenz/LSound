//                   ________________
// Tubo:      Vn ---|________________   -- aberto 
//

// Element size
lc = 0.05;

// Corners
Point(1) = { 0, 0,   0, lc};
Point(2) = { 1, 0,   0, lc};
Point(3) = { 1, 2*0.05,   0, lc};
Point(4) = { 0, 2*0.05,   0, lc};

// Edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Material
Physical Surface("Material,Ar,1,1.21,342.57,400.0") = {1};

// Boundary conditions - Open
Physical Curve("Open") = {2};

// Normal velocity (left)
Physical Curve("Vn,1E-3,100,0.0") = {4};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "tubo_teste.msh";
