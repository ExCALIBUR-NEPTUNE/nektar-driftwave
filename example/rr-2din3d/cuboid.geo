//=============================== Parameters ==================================
// Lengths and resolutions in each dimension
xsize = 100;
ysize = 100;
zsize = 36;
xmin=-xsize/2;
ymin=-ysize/2;
zmin=-zsize/2;
nx = 64;
ny = 64;
nz = 8; // Make dz ~= dx,dy
//=============================================================================

// Create a line of length <xsize> in the x-direction, with <nx> divisions
Point(1) = {xmin, ymin, zmin, 0.01};
Point(2) = {xmin+xsize, ymin, zmin, 0.01};
Line(1) = {1, 2};
Transfinite Line(1) = nx+1;

// Extrude split line into meshed square/rectangle
sq = Extrude {0,ysize,0} {Curve{1}; Layers{ny}; Recombine;};

// Extrude square/rectangle into a cuboid
cbd = Extrude {0,0,zsize} {Surface{sq[1]}; Layers{nz}; Recombine;};

// Define physical volume, surfaces for BCs
// Domain
Physical Volume(0) = {cbd[1]};
// Low-x side
Physical Surface(1) = {cbd[5]};
// High-x side
Physical Surface(2) = {cbd[3]};
// Low-y side
Physical Surface(3) = {cbd[2]};
// High-y side
Physical Surface(4) = {cbd[4]};
// Low-z side
Physical Surface(5) = {sq[1]};
// High-z side
Physical Surface(6) = {cbd[0]};