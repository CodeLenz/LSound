//
// Sala quadrada para estudarmos o BESO
//

// Element size
lc = 1e-1;

// Corners
Point(1) = { 0, 0,   0, lc};
Point(2) = { 1, 0,   0, lc};
Point(3) = { 1, 1,   0, lc};
Point(4) = { 0, 1,   0, lc};

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

// Aplica uma excitação de velocidade normal na parede da esquerda
// Physical Curve("Vn,1E-3,100.0,0.0") = {4};

// Aplica uma excitação de pressão imposta na parede da esquerda
// Só para testar a leitura de dados, ainda não implementado no 
// código
Physical Curve("Pressure,-1E-3,100.0,0.0") = {4};

// Nós para monitorar as pressões
// Vamos monitorar em TODOS os nós
Physical Surface("Probe") = {1};

// Nós para monitorar a SPL
Physical Curve("Target") = {2};

// Convert triangles to quads
Recombine Surface{:};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;

// Save the mesh
Save "sala_quadrada.msh";
