//=============================== Parameters ==================================
// Lengths and resolutions in each dimension
xsize = 2;
xmin = -1;
nx = 1000;
//=============================================================================

// Create a line in the x-direction of length <xsize>, with <nx> divisions
Point(1) = {xmin, 0, 0, 0.01};
Point(2) = {xmin+xsize, 0, 0, 0.01};
Line(0) = {1, 2};
Transfinite Line{0} = nx+1;

Physical Curve(0) = {0};
// Causes segfault in NekMesh; need to add manually
// Physical Point(1) = {1};
// Physical Point(2) = {2};
