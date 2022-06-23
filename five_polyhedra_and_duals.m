% Author: Osama M. Raisuddin

% This script creates 3D figures of the five regular polyhedra
% and figures with the duality demonstradted by overlaying them on each
% other




clc
clear 
close all


% Use the high quality renderer to ensure SVG plots are created
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')





% Choose the
% color,
% transparency,
% and viewing angle
% for the figures
color = [149, 121, 232]/255;
alpha = 0.8;
view_angle = [120,30];




%% The figures of the five platonic solids

% The simplex

% Write down the coordinates
tetrahedron_1 = [1, 1, 1;1, -1, -1;-1, 1, -1 ;-1, -1,  1];

figure
hold on
axis off
axis equal

% Generate the sets of vertices that form the faces of the solid
tetrahedron_vertex_sets = [[1,2,3];[1,2,4];[1,3,4];[2,3,4]];

% Generate faces of the solids
for i = 1:length(tetrahedron_vertex_sets)
    fill3(...
        tetrahedron_1(tetrahedron_vertex_sets(i,:),1),...
        tetrahedron_1(tetrahedron_vertex_sets(i,:),2),...
        tetrahedron_1(tetrahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)







% The dual of the first tetrahedron: another tetrahedron

% Write down the coordinates
tetrahedron_2 = [-1, -1, -1; -1, 1, 1; 1, -1, 1; 1,  1, -1];

figure
hold on
axis off
axis equal

% Generate faces of the solids
for i = 1:length(tetrahedron_vertex_sets)
    fill3(...
        tetrahedron_2(tetrahedron_vertex_sets(i,:),1),...
        tetrahedron_2(tetrahedron_vertex_sets(i,:),2),...
        tetrahedron_2(tetrahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)





% The octahedron

% Write down the coordinates
octahedron = [1,0,0;-1,0,0; 0,1,0;0,-1,0; 0,0,1;0,0,-1];

% Generate the sets of vertices that form the faces of the solid
octahedron_vertex_sets = [1,3,5; 1,3,6; 1,4,5; 1,4,6; 2,3,5; 2,3,6; 2,4,5; 2,4,6;];

figure
hold on
axis off
axis equal

% Generate faces of the solids
for i = 1:length(octahedron_vertex_sets)
    fill3(...
        octahedron(octahedron_vertex_sets(i,:),1),...
        octahedron(octahedron_vertex_sets(i,:),2),...
        octahedron(octahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)







% The cube

% Write down the coordinates
cube = [
     1, 1, 1;
     1, 1,-1;
     1,-1, 1;
     1,-1,-1;
    -1, 1, 1;
    -1, 1,-1;
    -1,-1, 1;
    -1,-1,-1;
     ];

% Generate the sets of vertices that form the faces of the solid
cube_vertex_sets = [
    1,2,4,3;
    1,5,7,3;
    1,2,6,5;
    8,7,5,6;
    8,4,3,7;
    8,4,2,6
    ];

figure
hold on
axis off
axis equal

% Generate faces of the solids
for i = 1:length(cube_vertex_sets)
    fill3(...
        cube(cube_vertex_sets(i,:),1),...
        cube(cube_vertex_sets(i,:),2),...
        cube(cube_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)








% The icosahedron

phi = (1+sqrt(5))/2;

% Write down the coordinates
icosahedron = [
    -phi, 1,   0;
     phi, 1,   0;
    -phi,-1,   0;
     phi,-1,   0;
       0, phi, 1;
       0,-phi, 1;
       0, phi,-1;
       0,-phi,-1;
      -1,   0, phi;
       1,   0, phi;
      -1,   0,-phi;
       1,   0,-phi
    ];

% Generate the connectivity matrix based on the distances between the
% vertices
icosahedron_connectivity_matrix = zeros(20,20);
for i = 1:12
    for j = i:12
        if sum((icosahedron(i,:)-icosahedron(j,:)).^2,2) == 4
            icosahedron_connectivity_matrix(i,j) = 1;
        end
    end
end

 % Generate the vertex sets using the connectivity matrix
icosahedron_vertex_sets = [];
for i = 1:12
    for j = i:12
        for k = j:12
            if icosahedron_connectivity_matrix(i,j) && icosahedron_connectivity_matrix(j,k) && icosahedron_connectivity_matrix(i,k)
                icosahedron_vertex_sets = [icosahedron_vertex_sets;[i,j,k]];
            end
        end
    end
end

figure
hold on
axis off
axis equal

% Generate faces of the solids
for i = 1:length(icosahedron_vertex_sets)
    fill3(...
        icosahedron(icosahedron_vertex_sets(i,:),1),...
        icosahedron(icosahedron_vertex_sets(i,:),2),...
        icosahedron(icosahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)






% The dodecahedron

phi = (1+sqrt(5))/2;

% Write down the coordinates
dodecahedron = [
     1, 1, 1;
     1, 1,-1;
     1,-1, 1;
     1,-1,-1;
    -1, 1, 1;
    -1, 1,-1;
    -1,-1, 1;
    -1,-1,-1;
    0, 1/phi, phi;
    0, 1/phi,-phi;
    0,-1/phi, phi;
    0,-1/phi,-phi;
     1/phi, phi,0;
     1/phi,-phi,0;
    -1/phi, phi,0;
    -1/phi,-phi,0;
     phi,0, 1/phi;
     phi,0,-1/phi;
    -phi,0, 1/phi;
    -phi,0,-1/phi;
    ];

% Generate the connectivity matrix based on the distances between the
% vertices
dodecahedron_connectivity_matrix = zeros(20,20);
for i = 1:20
    for j = 1:20
        if (sum((dodecahedron(i,:)-dodecahedron(j,:)).^2,2) < 2) && (sum((dodecahedron(i,:)-dodecahedron(j,:)).^2,2) ~= 0)
            dodecahedron_connectivity_matrix(i,j) = 1;
        end
    end
end


 % Generate the vertex sets using the connectivity matrix
dodecahedron_vertex_sets = [];
for i = 1:20
    for j = 1:20
        for k = 1:20
            for l = 1:20
                for m = 1:20
                    if dodecahedron_connectivity_matrix(i,j) && dodecahedron_connectivity_matrix(j,k) && dodecahedron_connectivity_matrix(k,l) && dodecahedron_connectivity_matrix(l,m) && dodecahedron_connectivity_matrix(m,i)
                        if isempty(dodecahedron_vertex_sets)
                            dodecahedron_vertex_sets = [dodecahedron_vertex_sets;[i,j,k,l,m]];
                        end
                        if ~ismember(perms([i,j,k,l,m]),dodecahedron_vertex_sets,'rows')
                            dodecahedron_vertex_sets = [dodecahedron_vertex_sets;[i,j,k,l,m]];
                        end
                    end
                end
            end
        end
    end
end

figure
hold on
axis off
axis equal

% Generate faces of the solids
for i = 1:length(dodecahedron_vertex_sets)
    fill3(...
        dodecahedron(dodecahedron_vertex_sets(i,:),1),...
        dodecahedron(dodecahedron_vertex_sets(i,:),2),...
        dodecahedron(dodecahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)














%% The figures of the duals

% The simplex and its dual (simplex)

figure
hold on
axis off
axis equal


% Generate the faces of the solid and its dual (different transparencies)
for i = 1:length(tetrahedron_vertex_sets)
    fill3(...
        tetrahedron_2(tetrahedron_vertex_sets(i,:),1),...
        tetrahedron_2(tetrahedron_vertex_sets(i,:),2),...
        tetrahedron_2(tetrahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',0)
end
for i = 1:length(tetrahedron_vertex_sets)
    fill3(...
        tetrahedron_1(tetrahedron_vertex_sets(i,:),1),...
        tetrahedron_1(tetrahedron_vertex_sets(i,:),2),...
        tetrahedron_1(tetrahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)








% The dodecahedron and dual

figure
hold on
axis off
axis equal

% Generate the faces of the solid and its dual (different transparencies)
for i = 1:length(octahedron_vertex_sets)
    fill3(...
        octahedron(octahedron_vertex_sets(i,:),1),...
        octahedron(octahedron_vertex_sets(i,:),2),...
        octahedron(octahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end
for i = 1:length(cube_vertex_sets)
    fill3(...
        cube(cube_vertex_sets(i,:),1),...
        cube(cube_vertex_sets(i,:),2),...
        cube(cube_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',0)
end

% Set the viewing angle
view(view_angle)







% The icosahedron and its dual (dodecahedron)

figure
hold on
axis off
axis equal

% Generate the faces of the solid and its dual (different transparencies)
for i = 1:length(icosahedron_vertex_sets)
    fill3(...
        icosahedron(icosahedron_vertex_sets(i,:),1),...
        icosahedron(icosahedron_vertex_sets(i,:),2),...
        icosahedron(icosahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',0)
end
for i = 1:length(dodecahedron_vertex_sets)
    fill3(...
        dodecahedron(dodecahedron_vertex_sets(i,:),1),...
        dodecahedron(dodecahedron_vertex_sets(i,:),2),...
        dodecahedron(dodecahedron_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)

