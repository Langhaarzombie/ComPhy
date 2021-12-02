cl1 = 10.0;

// Box points

// COIL //
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


// EXTRUSIONS
coil[] = Extrude{0, 0, 100} {
    Surface{1};
};
