cl1 = 1;
Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {-1, 0, 0, cl1};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Line Loop(4) = {1, 2};
Plane Surface(4) = {4};
