// COIL //
cl1 = 10.0;
// Points
// Inner
Point(1) = {35, 35, -50, cl1};
Point(2) = {35, -35, -50, cl1};
Point(3) = {-35, -35, -50, cl1};
Point(4) = {-35, 35, -50, cl1};
Point(5) = {-5, -35, -50, cl1}; // used for the gap
Point(6) = {5, -35, -50, cl1}; // used for the gap

// Outer
Point(7) = {45, 45, -50, cl1};
Point(8) = {45, -45, -50, cl1};
Point(9) = {-45, -45, -50, cl1};
Point(10) = {-45, 45, -50, cl1};
Point(11) = {-5, -45, -50, cl1}; // used for the gap
Point(12) = {5, -45, -50, cl1}; // used for the gap

// Lines
// Inner
Line(1) = {1, 2};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {2, 6};
Line(5) = {3, 5};

// Outer
Line(6) = {7, 8};
Line(7) = {7, 10};
Line(8) = {10, 9};
Line(9) = {8, 12};
Line(10) = {9, 11};

// Inner - Outer
Line(11) = {5, 11};
Line(12) = {6, 12};

Line Loop(1) = {1, 4, 12, -9, -6, 7, 8, 10, -11, -5, -3, -2};

// Surface
Plane Surface(1) = {1};

// BOX //
cl2 = 5.0;
// Points
Point(13) = {150, 150, -150, cl2};
Point(14) = {150, -150, -150, cl2};
Point(15) = {-150, -150, -150, cl2};
Point(16) = {-150, 150, -150, cl2};

// Lines
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};

Line Loop(2) = {13, 14, 15, 16};

// Surface
Plane Surface(2) = {2};

// EXTRUSIONS
coil[] = Extrude{0, 0, 100} {
    Surface{1};
};

cube[] = Extrude{0, 0, 300} {
    Surface{2};
};

// Delete unused volumes (else there is an error with volumes overlappting)
Delete{
    Volume{coil[1]};
    Volume{cube[1]};
}

Surface Loop(1) = {
    1, coil[0], coil[2], coil[3], coil[4], coil[5], coil[6], coil[7], coil[8], coil[9], coil[10], coil[11], coil[12], coil[13]
};

Surface Loop(2) = {
    2, cube[0], cube[2], cube[3], cube[4], cube[5]
};

// VOLUMES
Volume(100) = {1};
Volume(200) = {2, 1}; // first argument is the exterior, second the holes

Physical Volume(1) = {100};
Physical Volume(2) = {200};
