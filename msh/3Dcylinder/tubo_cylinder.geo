//
// 3D tube
//

SetFactory("OpenCASCADE");

// Tamanho dos elementos
//Mesh.CharacteristicLength{:} = 0.25;

// Volume (box 3D)
Cylinder(1) = {0, 0, 0, 1, 0, 0, 0.02, 2*Pi};

// Material
Physical Volume("Material,Ar,1,1.028,340.0,400.0") = {1};


Mesh 2;
RefineMesh;
Mesh 3;
