// Gmsh project created on Mon Apr 30 11:43:39 2018
SetFactory("OpenCASCADE");
//+
air = 1;
capH = 0.2;
capW = 0.2;

Rectangle(1) = {-air/2, -air/2, 0, air, air, 0};
Rectangle(2) = {-capW/2, -capH/2, 0, capW/2, capH, 0};
Rectangle(3) = {0, -capH/2, 0, capW/2, capH, 0};

BooleanDifference{Surface{1};Delete;}{Surface{2};Surface{3};}
Coherence;
Characteristic Length{8,17,11,10,16,5} = 0.01;
Characteristic Length{12,13,14,15} = 0.05;

Physical Surface("M:Left") = {2};
Physical Surface("M:Right") = {3};
Physical Surface("M:Air") = {1};

Physical Line("BD:Left_Plate") = {8};
Physical Line("BD:Right_Plate") = {10};

Mesh 2;

Mesh.Format = 1;

Save "planar_gmsh.msh";
