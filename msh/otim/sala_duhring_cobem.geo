//
// Sala do artigo da Maria Duhring
//

SetFactory("OpenCASCADE");

// Element size
lc = 0.2;

// Tamanho de elemento 
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// Pontos do domínio
Point(1) = { 0  ,  0,   0, lc};
Point(2) = { 18 ,  0,   0, lc};
Point(3) = { 18 ,  8,   0, lc};
Point(4) = { 0  ,  8,   0, lc};
Point(5) = { 18 ,  9,   0, lc};
Point(6) = { 0  ,  9,   0, lc};

// Linhas da parte que não é de projeto
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Linhas da parte de projeto
Line(5) = {3, 5};
Line(6) = {5, 6};
Line(7) = {6, 4};

// Região de excitação 
Circle(8) = {2, 2, 0, 0.1, 0, 2*Pi};

// Região target
Circle(9) = {16, 2, 0, 1, 0, 2*Pi};

// Loop para a área que não é de projeto
Curve Loop(10) = {1, 2, 3, 4}; 

// Loop para a área de projeto
Curve Loop(20) = {-3, 5, 6, 7};

// Curve loops para os círculos
Curve Loop(30) = {8};
Curve Loop(40) = {9};

// Superfície que não é de projeto
Plane Surface(100) = {10};
Curve{8, 9} In Surface{100};

// Superfície de projeto
Plane Surface(200) = {20};

// Adicionando os circulos na superficie que não é de projeto
Curve{8, 9} In Surface{100};

// Material da sala (Ar)
Physical Surface("Material,Ar,1,1.204,343.328,400.0") = {100};

// Material do teto (Aluminio)
Physical Surface("Material,Al,2,2643.0,5098.3516,400.0") = {200};

// Aplica uma excitação de velocidade normal na parede da esquerda
Physical Curve("Vn,1E-3,34.56,0.0") = {8};

// A região 1 (Ar) não é de projeto
Physical Surface("Fixed,0.0") = {100};

// Aplica uma excitação de pressão imposta na parede da esquerda
// Só para testar a leitura de dados, ainda não implementado no 
// código
//Physical Curve("Pressure,1E1,100.0,0.0") = {4};

// Nós para monitorar as pressões
// Vamos monitorar em TODOS os nós
// Physical Surface("Probe") = {1};

// Nós para monitorar a SPL
Physical Curve("Target") = {9};

// Convert triangles to quads
Recombine Surface{200};

// Better quad algorithm
Mesh.Algorithm = 8;

// Build mesh
Mesh 2;
Save "sala_duhring_cobem.msh";
