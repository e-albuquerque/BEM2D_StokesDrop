// dad_1.geo - Straight plane microchannel (L = 10H, H = 1) with 1 drop (a/H = 0.50)
lc = 0.15;

// Outer channel contour vertices (counter-clockwise)
Point(1) = {0.0,  0.0, 0.0, lc};
Point(2) = {10.0, 0.0, 0.0, lc};
Point(3) = {10.0, 1.0, 0.0, lc};
Point(4) = {0.0,  1.0, 0.0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Drop contour: center (0.5, 0.5), radius R = 0.25 (clockwise starting at theta = 0)
Point(5) = {0.5,  0.5,   0.0, lc};
Point(6) = {0.75, 0.5,   0.0, lc};
Point(7) = {0.375, 0.5 - 0.21650635094610965, 0.0, lc};
Point(8) = {0.375, 0.5 + 0.21650635094610965, 0.0, lc};

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 6};

// Boundary discretization (470 quadratic elements: 440 on channel + 30 on drop)
Transfinite Curve {1, 3} = 201;
Transfinite Curve {2, 4} = 21;
Transfinite Curve {5, 6, 7} = 11;

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7};
Plane Surface(1) = {1, 2};
Plane Surface(2) = {-2};

Physical Curve("wall", 1)   = {1, 3};
Physical Curve("outlet", 2) = {2};
Physical Curve("inlet", 3)  = {4};
Physical Curve("drop_1", 4) = {5, 6, 7};
Physical Surface("fluid", 10) = {1};
Physical Surface("drop_surf_1", 11) = {2};
