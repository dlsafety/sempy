l = 1;
r2 = 0.57735;
n = 5;

Point(1) = {0,0,0,l};
Point(102) = {r2,r2,-r2,l};
Point(103) = {-r2,r2,-r2,l};
Point(104) = {-r2,-r2,-r2,l};
Point(105) = {r2,-r2,-r2,l};
Circle(29) = {103,1,102};
Circle(30) = {102,1,105};
Circle(31) = {105,1,104};
Circle(32) = {104,1,103};
Line Loop(33) = {29,30,31,32};
Ruled Surface(34) = {33};
Rotate { {1,0,0},{0,0,0}, Pi/2 } { Duplicata{ Surface{34}; } }
Rotate { {1,0,0},{0,0,0}, Pi } { Duplicata{ Surface{34}; } }
Rotate { {1,0,0},{0,0,0}, 3*Pi/2 } { Duplicata{ Surface{34}; } }
Rotate { {0,1,0},{0,0,0}, Pi/2 } { Duplicata { Surface{34}; } }
Rotate { {0,1,0},{0,0,0}, -Pi/2 } { Duplicata { Surface{34}; } }

Surface Loop(101) = {35, 40, 45, 34, 51, 50};
Volume(102) = {101};

Transfinite Line {42, 44, 30, 32, 36, 41, 29, 31, 39, 49, 37, 47} = n;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";


