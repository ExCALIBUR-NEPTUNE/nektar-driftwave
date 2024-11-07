D=100;
r=D/2;
// Choose size to generate approx 128^2 elements
element_size=1.07;

Point(5) = { 0,  0, 0, 0.3};
Point(6) = { r,  0, 0, 0.3};
Point(7) = { 0,  r, 0, 0.3};
Point(8) = {-r,  0, 0, 0.3};
Point(9) = { 0, -r, 0, 0.3};

// Draw circle in four segments
Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

// Create line loop and surface
Line Loop(1) = {5,6,7,8};
Plane Surface(1) = {1};

// Label boundary
Physical Line(1) = {5,6,7,8};
// Label domain
Physical Surface(0) = {1};

Mesh.CharacteristicLengthMin = element_size;