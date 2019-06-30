//SetFactory("OpenCASCADE");
// lc = 0.009;
// lc = 0.005;
lc = 0.05;
lc_in = 0.05;
lc_out = 0.05;
// lc = 0.01;
// lc = 0.007;

r = 0.4;
h1 = 0.1;
l = 4;

Point(1) = {0, -r, 0, lc_in};
Point(2) = {0, r, 0, lc_in};
Line(1) = {1,2};
Point(3) = {l/4, r, 0 , lc};
Point(4) = {l/8+l/4, r+h1/2, 0 , lc};
Point(5) = {l/2, r + h1, 0, lc};

BSpline(2) = {2, 3, 4, 5};

h2 = 0.5;

Point(6) = {l/2+l/8, r + 3*h1/2, 0, lc};

Point(7) = {13/16*l, r + 3*h1/2+h2/2,0, lc};
Point(8) = {l, r + 3*h1/2+h2,0, lc_out};
BSpline(3) = {5, 6, 7, 8};

nnorm = Sqrt( (3/16*l)*(3/16*l) + (h2/2)*(h2/2));
n1 = (3/16*l)/nnorm;
n2 = (h2/2)/nnorm;

r1 = 0.5;

t1 = n2;
t2 = -n1;

Point(9) = {l + t1*r1, r + 3*h1/2+h2 + t2*r1,0, lc_out};
//+
Line(4) = {8, 9};

Point(10) = {13/16*l + t1*r1, r + 3*h1/2+h2/2 + t2*r1,0, lc};
Point(11) = {l/2+l/8 + t1*r1, r + 3*h1/2 + t2*r1, 0, lc};

h3 = 0.00;
Point(13) = {l/4, -r, 0 , lc};
Point(14) = {l/8+l/4, -(r+h3/2), 0 , lc};
Point(15) = {l/2, -(r + h3), 0, lc};

BSpline(5) = {1, 13, 14, 15};

Point(16) = {l/2+l/8, -(r + 3*h3/2), 0, lc};

h4 = 0.1;
Point(17) = {13/16*l, -(r + 3*h3/2+h4/2),0, lc};
Point(18) = {l, -(r + 3*h3/2+h4),0, lc_out};
BSpline(6) = {15, 16, 17, 18};

nnorm = Sqrt( (3/16*l)*(3/16*l) + (h4/2)*(h4/2));
n1 = (3/16*l)/nnorm;
n2 = -(h4/2)/nnorm;

r2 = 0.3;

t1_ = -n2;
t2_ = +n1;

Point(19) = {l + t1_*r2, -(r + 3*h3/2+h4) + t2_*r2,0, lc_out};
//+
Line(7) = {18, 19};

Point(20) = {13/16*l+ t1_*r2, -(r + 3*h3/2+h4/2)+ t2_*r2,0, lc};
Point(21) = {l/2+l/8+ t1_*r2, -(r + 3*h3/2)+ t2_*r2, 0, lc};


Point(12) = {(l/2+l/8 + t1*r1 + l/2+l/8+ t1_*r2)/2, (r + 3*h1/2 + t2*r1 -(r + 3*h3/2)+ t2_*r2)/2, 0, lc};
//+
BSpline(8) = {19, 20, 21, 12};
//+
BSpline(9) = {12, 11, 10, 9};
//+
Line Loop(1) = {6, 7, 8, 9, -4, -3, -2, -1, 5};
//+
Plane Surface(1) = {1};
//+
Physical Line("Wall_down", 1) = {5, 6};
//+
Physical Line("Outlet_1", 2) = {7};
//+
Physical Line("Wall_right", 3) = {8, 9};
//+
Physical Line("Outlet_2", 4) = {4};
//+
Physical Line("Wall_up", 5) = {3, 2};
//+
Physical Line("Inlet", 6) = {1};
//+
Physical Surface("S", 7) = {1};
