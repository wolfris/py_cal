cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {0.7, 0, 0, 1};
Point(3) = {1, 0, 0, 1};
Point(4) = {0, 1, 0, 1};
Point(5) = {0, -1, 0, 1};
Point(6) = {-0.5, 0, 0, 1};
Point(7) = {-0.4, 0, 0, 1};
Point(8) = {-0.6, 0, 0, 1};
Point(9) = {0.2, -0.3, 0, 1};
Point(10) = {0.4, -0.3, 0, 1};
Point(11) = {0.6, -0.3, 0, 1};
Line(1) = {2, 1};
Line(2) = {2, 3};
Line(3) = {1, 4};
Circle(4) = {3, 1, 5};
Circle(5) = {4, 1, 5};
Circle(9) = {7, 6, 8};
Circle(10) = {8, 6, 7};
Circle(11) = {11, 10, 9};
Circle(12) = {9, 10, 11};
Line Loop(207) = {3, 5, -4, -2, 1, -10, -9, -11, -12};
Plane Surface(207) = {207};
Physical Line(200) = {3};
Physical Line(201) = {1, 2, 4, 5, 9, 10, 11, 12};
Physical Surface(208) = {207};
