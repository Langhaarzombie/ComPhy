cl = 0.2;
lz = 2;

Point(1) = {0,  0, 0, cl};
Point(2) = {5,  1, 0, cl};
Point(3) = {10, 0, 0, cl};
Point(4) = {0,  4, 0, cl};
Point(5) = {5,  3, 0, cl};
Point(6) = {10, 4, 0, cl};

Line(1) = {1,4};
Line(2) = {3,6};
Spline(3) = {1,2,3};
Spline(4) = {4,5,6};

Line Loop(1) = {1,4,-2,-3};
Plane Surface(1) = {1};

s0[] = Extrude{0, 0, lz} {Surface{1};};
