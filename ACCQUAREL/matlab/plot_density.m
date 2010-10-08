% import density.
bigmat = dlmread('../density');
x = bigmat(:, 1);
y = bigmat(:, 2);
z = bigmat(:, 3);
d = bigmat(:, 4);
N = round(numel(x)^(1/3)); % number of points in each dimension (assume same mesh in the three dimensions)
L = max(x); % mesh is a cube [-L, L]^3
M = max(abs(d)); % maximum value of density

x = linspace(-L, L, N);
y = x;
z = x;
[X, Y, Z] = meshgrid(x, y, z); % generates the internal mesh

D = reshape(d, N, N, N); % D is a 3D matrix
clf
%slice(X, Y, Z, (D), 0, 0, 0); % to plot slices of the volume
ms = [M/1.1 M/2 M/5 M/10 M/20 M/50 M/100 M/200 M/400]; % which isosurfaces to plot
hs = arrayfun(@(m) patch(isosurface(X, Y, Z, D, m, log(D))), ms); % plot each isosurface
isonormals(X,Y,Z,D,D);
shading flat % shading faceted for the mesh
set(hs, 'FaceAlpha', .4);
axis equal
axis off

xlabel 'x'
ylabel 'y'
zlabel 'z'

view(3)
camlight 'right'
