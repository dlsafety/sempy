// Gmsh project created on Sun Feb 12 17:29:38 2017
Point(1) = {-1, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {0,-1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, 0, 0, 1.0};

Point(6) = {-.5, 0, 0, 1.0};
Point(7) = { .5, 0, 0, 1.0};
Point(8) = {0, -.5, 0, 1.0};
Point(9) = {0,  .5, 0, 1.0};

Circle(1) = {4, 5, 2};
Circle(2) = {2, 5, 3};
Circle(3) = {3, 5, 1};
Circle(4) = {1, 5, 4};
Line(5) = {9, 7};
Line(6) = {7, 8};
Line(7) = {8, 6};
Line(8) = {6, 9};
Line(9) = {4, 9};
Line(10) = {2, 7};
Line(11) = {3, 8};
Line(12) = {1, 6};
Line Loop(13) = {4, 9, -8, -12};
Plane Surface(14) = {13};
Line Loop(15) = {12, -7, -11, 3};
Plane Surface(16) = {15};
Line Loop(17) = {11, -6, -10, 2};
Plane Surface(18) = {17};
Line Loop(19) = {10, -5, -9, 1};
Plane Surface(20) = {19};
Line Loop(21) = {8, 5, 6, 7};
Plane Surface(22) = {21};

Extrude {0, 0, 2} {
  Line{4, 9, 1, 5, 10, 6, 2, 7, 11, 3, 12, 8};
}
Line Loop(71) = {23, 27, -67, -63};
Plane Surface(72) = {71};
Line Loop(73) = {59, 63, -51, -55};
Plane Surface(74) = {73};
Line Loop(75) = {55, -43, -39, 47};
Plane Surface(76) = {75};
Line Loop(77) = {35, -39, -31, 27};
Plane Surface(78) = {77};
Line Loop(79) = {35, 43, 51, 67};
Plane Surface(80) = {79};

Surface Loop(81) = {62, 16, 74, 66, 54, 58};
Volume(82) = {81};
Surface Loop(83) = {58, 50, 18, 76, 46, 42};
Volume(84) = {83};
Surface Loop(85) = {42, 38, 20, 34, 78, 30};
Volume(86) = {85};
Surface Loop(87) = {72, 26, 14, 70, 30, 66};
Volume(88) = {87};
Surface Loop(89) = {54, 38, 22, 80, 70, 46};
Volume(90) = {89};

Physical Surface(91) = {72, 78, 34, 20, 18, 50, 74, 76, 80, 62, 26, 14, 16, 22};
Physical Volume(92) = {84, 90, 88, 86, 82};

// Transfinite Line "*";
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
