//
//
// Arquivo utilizado para testar a resposta de um tubo fechado - fechado
// submetido a uma excitação de velocidade normal na face fechada da esquerda.
// O material é ar, com as propriedades rho = 1.21 e c = 342.57 (k = 1.42E5)
//

//                   ________________
// Tubo:      Vn ---|________________| 
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

// Nós para monitorar a SPL
Physical Curve("Target") = {2};

// Normal velocity (left)
Physical Curve("Vn,1E-3,100,0.0") = {4};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "tubo_fechado_fechado.msh";
