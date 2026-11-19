// dad_2.geo - Constricted microcapillary channel (L = 3H, H/h = 2) with 1 drop (a/H = 0.60)
lc = 0.1;

// Outer channel contour vertices (counter-clockwise)
Point(1)  = {0.0, 0.0,  0.0, lc};
Point(2)  = {1.0, 0.0,  0.0, lc};
Point(3)  = {1.0, 0.25, 0.0, lc};
Point(4)  = {2.0, 0.25, 0.0, lc};
Point(5)  = {2.0, 0.0,  0.0, lc};
Point(6)  = {3.0, 0.0,  0.0, lc};
Point(7)  = {3.0, 1.0,  0.0, lc};
Point(8)  = {2.0, 1.0,  0.0, lc};
Point(9)  = {2.0, 0.75, 0.0, lc};
Point(10) = {1.0, 0.75, 0.0, lc};
Point(11) = {1.0, 1.0,  0.0, lc};
Point(12) = {0.0, 1.0,  0.0, lc};

// Channel boundary segments
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 5};
Line(5)  = {5, 6};
Line(6)  = {6, 7};
Line(7)  = {7, 8};
Line(8)  = {8, 9};
Line(9)  = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

// Drop contour: center (0.5, 0.5), radius R = 0.3 (clockwise starting at theta = 0)
Point(13) = {0.5, 0.5, 0.0, lc};
Point(14) = {0.8, 0.5, 0.0, lc};
Point(15) = {0.5, 0.2, 0.0, lc};
Point(16) = {0.2, 0.5, 0.0, lc};
Point(17) = {0.5, 0.8, 0.0, lc};

Circle(13) = {14, 13, 15};
Circle(14) = {15, 13, 16};
Circle(15) = {16, 13, 17};
Circle(16) = {17, 13, 14};

// Boundary discretization (220 quadratic elements: 180 on channel + 40 on drop)
Transfinite Curve {1, 3, 5, 7, 9, 11} = 21;
Transfinite Curve {2, 4, 8, 10} = 6;
Transfinite Curve {6, 12} = 21;
Transfinite Curve {13, 14, 15, 16} = 11;

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Curve Loop(2) = {13, 14, 15, 16};
Plane Surface(1) = {1, 2};
Plane Surface(2) = {-2};

Physical Curve("wall", 1)   = {1, 2, 3, 4, 5, 7, 8, 9, 10, 11};
Physical Curve("outlet", 2) = {6};
Physical Curve("inlet", 3)  = {12};
Physical Curve("drop_1", 4) = {13, 14, 15, 16};
Physical Surface("fluid", 10) = {1};
Physical Surface("drop_surf_1", 11) = {2};
