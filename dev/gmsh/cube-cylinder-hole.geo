// Gmsh project created on Sat Feb 11 08:03:57 2017

Point(1) = {0, -0.5, 0, 1.0};
Point(2) = {0, 0.5, 0, 1.0};
Point(3) = {0.5, 0, 0, 1.0};
Point(4) = {-0.5, 0, 0, 1.0};
Point(5) = {-0.1, 0, 0, 1.0};
Point(6) = {0.1, 0, 0, 1.0};
Point(7) = {0, 0.1, 0, 1.0};
Point(8) = {0, -0.1, 0, 1.0};
Line(1) = {4, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {5, 4};
Line(6) = {8, 8};
Line(7) = {8, 1};
Line(8) = {6, 3};
Line(9) = {7, 2};
Point(9) = {0, 0, 0, 1.0};
Circle(10) = {7, 9, 5};
Circle(11) = {7, 9, 6};
Circle(12) = {6, 9, 8};
Circle(13) = {8, 9, 5};
Line Loop(14) = {1, -9, 10, 5};
Plane Surface(15) = {14};
Line Loop(16) = {5, -4, -7, 13};
Plane Surface(17) = {16};
Line Loop(18) = {12, 7, -3, -8};
Plane Surface(19) = {18};
Line Loop(20) = {9, 2, -8, -11};
Plane Surface(21) = {20};


Extrude {0, 0, 1} {
  Line{1, 5, 4, 7, 3, 8, 2, 9, 10, 13, 12, 11};
}
Line Loop(70) = {30, -26, -58, 34};
Plane Surface(71) = {70};
Line Loop(72) = {34, -38, -42, 62};
Plane Surface(73) = {72};
Line Loop(74) = {66, 42, -46, -50};
Plane Surface(75) = {74};
Line Loop(76) = {50, -22, -26, -54};
Plane Surface(77) = {76};


Surface Loop(78) = {41, 19, 65, 73, 37, 45};
Volume(79) = {78};
Surface Loop(80) = {61, 17, 33, 71, 37, 29};
Volume(81) = {80};
Surface Loop(82) = {25, 15, 57, 77, 29, 53};
Volume(83) = {82};
Surface Loop(84) = {69, 21, 49, 75, 45, 53};
Volume(85) = {84};


// Transfinite Line "*";
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
