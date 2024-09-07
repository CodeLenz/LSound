//
// 3D tube
//

SetFactory("OpenCASCADE");

// Tamanho dos elementos
//Mesh.CharacteristicLength{:} = 0.25;

// Volume (box 3D)
Box(1) = {0, 0, 0, 1, 0.1, 0.1};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

// Material
Physical Volume("Material,Ar,1,1.028,340.0,400.0") = {1};

// Normal velocity (left)
Physical Surface("Vn,1E-3,340.0,0.0") = {1};

// Probe
Physical Curve("Probe") = {10};

// Atenuation
Physical Surface("Yn,1.0") = {1,2};

// Better quad algorithm
Mesh.Algorithm = 8;

// Para gerar a malha
Mesh.Algorithm3D = 9;

// Build mesh
Mesh 3;

// Save the mesh
Save "tubo_3d.msh";
