//
// Tube 
//

// Element size
lc = 0.001;

// Variables

// Tube
L = 0.1;
r = 0.01;

// Neck
Lp = 0.02;
rp = 0.01;

// Cavity
rc = 3E-2;
Hc = 5E-2;


//  Bottom line
Point(1) = { 0,        0,     0, lc};
Point(2) = { L,        0,     0, lc};
Point(3) = { L+2*rp,   0,     0, lc};
Point(4) = { L+2*rp+4*L, 0,     0, lc};

// top of the Tube
Point(5) = { 0,        2*r,     0, lc};
Point(6) = { L,        2*r,     0, lc};
Point(7) = { L+2*rp,   2*r,     0, lc};
Point(8) = { L+2*rp+4*L, 2*r,     0, lc};

// Bottom of Cavity
Point(9)  = { L-rc,        2*r+Lp,     0, lc};
Point(10) = { L,           2*r+Lp,     0, lc};
Point(11) = { L+2*rp,      2*r+Lp,     0, lc};
Point(12) = { L+2*rp+rc,   2*r+Lp,     0, lc};

// Top of Cavity
Point(13)  = { L-rc,        2*r+Lp+Hc,     0, lc};
Point(14)  = { L,           2*r+Lp+Hc,     0, lc};
Point(15)  = { L+2*rp,      2*r+Lp+Hc,     0, lc};
Point(16)  = { L+2*rp+rc,   2*r+Lp+Hc,     0, lc};


// Edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Line(4) = {4, 8};
Line(5) = {8, 7};
Line(6) = {7, 11};
Line(7) = {11, 12};
Line(8) = {12, 16};
Line(9) = {16, 15};
Line(10) = {15, 14};
Line(11) = {14, 13};
Line(12) = {13, 9};
Line(13) = {9, 10};
Line(14) = {10, 6};
Line(15) = {6, 5};
Line(16) = {5, 1};


// Surface
Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
Plane Surface(1) = {1};


// Material
Physical Surface("Material,Ar,1,1.028,340.0,400.0") = {1};

// Normal velocity (left)
Physical Curve("Vn,1E-3,741.0,0.0") = {16};

// Probe nodes
//Physical Curve("Probe") = {2};

//Physical Curve("Open") = {4};

// Atenuation
Physical Curve("Yn,1.0") = {4};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

