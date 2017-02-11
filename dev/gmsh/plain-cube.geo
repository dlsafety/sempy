
// Gmsh project created on Sat Feb 11 13:01:09 2017

Point(1) = {-1, -1, -1, 1.0};
Point(2) = { 1, -1, -1, 1.0};
Point(3) = {-1,  1, -1, 1.0};
Point(4) = { 1,  1, -1, 1.0};
Point(5) = {-1, -1,  1, 1.0};
Point(6) = { 1, -1,  1, 1.0};
Point(7) = {-1,  1,  1, 1.0};
Point(8) = { 1,  1,  1, 1.0};


Line(1) = {6, 8};
Line(2) = {8, 4};
Line(3) = {4, 2};
Line(4) = {2, 6};
Line(5) = {6, 5};
Line(6) = {5, 7};
Line(7) = {7, 3};
Line(8) = {3, 1};
Line(9) = {1, 5};
Line(10) = {7, 8};
Line(11) = {4, 3};
Line(12) = {2, 1};
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {10, -1, 5, 6};
Plane Surface(16) = {15};
Line Loop(17) = {4, 5, -9, -12};
Plane Surface(18) = {17};
Line Loop(19) = {3, 12, -8, -11};
Plane Surface(20) = {19};
Line Loop(21) = {2, 11, -7, 10};
Plane Surface(22) = {21};
Line Loop(23) = {8, 9, 6, 7};
Plane Surface(24) = {23};
Surface Loop(25) = {22, 14, 16, 18, 24, 20};
Volume(26) = {25};


// Transfinite Line "*";
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
