// Gmsh project created on Sun Feb 12 19:41:00 2017
Point(1) = {-1, -1, 0, 1.0};
Point(2) = {1, -1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {-1, 1, 0, 1.0};
Point(5) = {0, 0, 0, 1.0};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Circle(4) = {4, 5, 3};
Line Loop(5) = {1, 2, 3, -4};
Plane Surface(6) = {5};

Extrude {0, 0, 2} {
  Line{1, 4, 3, 2};
}

Line Loop(23) = {11, -15, -19, -7};
Plane Surface(24) = {23};

Surface Loop(25) = {24, 14, 6, 10, 22, 18};
Volume(26) = {25};

// Transfinite Line "*";
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
