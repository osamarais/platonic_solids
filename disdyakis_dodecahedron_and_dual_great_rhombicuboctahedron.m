% Author: Osama M. Raisuddin

% This script creates 3D figures of the disdyakis dodecahedron and its
% dual, the great rhombicuboctahedron



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







% Disdyakis dodecahedron

% Write down the coordinates
a = 1/(1+2*sqrt(2));
b = 1/(2+3*sqrt(2));
c = 1/sqrt(27+18*sqrt(2));
disd = [
     a, 0, 0;
    -a, 0, 0;
     0, a, 0;
     0,-a, 0;
     0, 0, a;
     0, 0,-a;
    
     b, b, 0;
     b,-b, 0;
    -b, b, 0;
    -b,-b, 0;
    
     0, b, b;
     0, b,-b;
     0,-b, b;
     0,-b,-b;
    
     b, 0, b;
     b, 0,-b;
    -b, 0, b;
    -b, 0,-b;
    
    
     c, c, c;
     c, c,-c;
     c,-c, c;
     c,-c,-c;
     
    -c, c, c;
    -c, c,-c;
    -c,-c, c;
    -c,-c,-c;
     
];

% Normalize the lengths
for i = 1:length(disd)
    disd(i,:) = disd(i,:)/norm(disd(i,:));
end

% Generate the connectivity matrix based on the distances between the
% vertices
disd_conn_mat = [];
for i = 1:length(disd)
    for j = 1:length(disd)
        if i==j
            continue
        end
        if (sum((disd(i,:)-disd(j,:)).^2,2) < 0.8454) && (sum((disd(i,:)-disd(j,:)).^2,2) ~= 0)
            disd_conn_mat(i,j) = 1;
        end
    end
end

% Generate the vertex sets using the connectivity matrix
disd_vert_sets = [];
for i = 1:length(disd_conn_mat)
    for j = 1:length(disd_conn_mat)
        for k = 1:length(disd_conn_mat)
                

                
                    if disd_conn_mat(i,j) & disd_conn_mat(j,k) & disd_conn_mat(k,i)
                        if isempty(disd_vert_sets)
                            disd_vert_sets = [disd_vert_sets;[i,j,k]];
                        end
                        if ~ismember(perms([i,j,k]),disd_vert_sets,'rows')
                            disd_vert_sets = [disd_vert_sets;[i,j,k]];
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
for i = 1:length(disd_vert_sets)
    fill3(...
        disd(disd_vert_sets(i,:),1),...
        disd(disd_vert_sets(i,:),2),...
        disd(disd_vert_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)





% Great rhombicuboctahedron

% Generate the coordinates

% Write down the base, of whose permutations the vertices are formed
base = [
     1, (1+sqrt(2)), (1+2*sqrt(2));
     1, (1+sqrt(2)),-(1+2*sqrt(2));
     1,-(1+sqrt(2)), (1+2*sqrt(2));
     1,-(1+sqrt(2)),-(1+2*sqrt(2));
    -1, (1+sqrt(2)), (1+2*sqrt(2));
    -1, (1+sqrt(2)),-(1+2*sqrt(2));
    -1,-(1+sqrt(2)), (1+2*sqrt(2));
    -1,-(1+sqrt(2)),-(1+2*sqrt(2));
    ];

% Make a single entry to avoid ismember() from failing
cubo = [1, (1+sqrt(2)), (1+2*sqrt(2));];
% Generate the remaining permutations
for i = 1:length(base)
    myperms = perms(base(i,:));
    for j = 1:length(myperms)
        if ~ismember(myperms(j,:),cubo,'rows')
            cubo = [cubo; myperms(j,:)];
        end
    end
end

% Generate the connectivity matrix based on the distances between the
% vertices
cubo_conn_mat = zeros(length(cubo),length(cubo));
for i = 1:length(cubo)
    for j = 1:length(cubo)
        if (sum((cubo(i,:)-cubo(j,:)).^2,2) <= 4.01) && (sum((cubo(i,:)-cubo(j,:)).^2,2) ~= 0)
            cubo_conn_mat(i,j) = 1;
        end
    end
end

% Normalize the coordinates
for i = 1:length(cubo)
    cubo(i,:) = cubo(i,:)/norm(cubo(i,:));
end


% Generate the vertex sets using the connectivity matrix
cubo_vertex_sets_square = [];
cubo_vertex_sets_hex = [];
cubo_vertex_sets_oct = [];
% Squares
for i = 1:length(cubo_conn_mat)
    for j = 1:length(cubo_conn_mat)
        for k = 1:length(cubo_conn_mat)
            for l = 1:length(cubo_conn_mat)
                
                if i==k | j==l
                    continue
                end
                
                    if cubo_conn_mat(i,j) & cubo_conn_mat(j,k) & cubo_conn_mat(k,l) & cubo_conn_mat(l,i)
                        if isempty(cubo_vertex_sets_square)
                            cubo_vertex_sets_square = [cubo_vertex_sets_square;[i,j,k,l]];
                        end
                        if ~ismember(perms([i,j,k,l]),cubo_vertex_sets_square,'rows')
                            cubo_vertex_sets_square = [cubo_vertex_sets_square;[i,j,k,l]];
                        end
                end
            end
        end
    end
end
% Hexagons
for i = 1:length(cubo_conn_mat)
    for j = 1:length(cubo_conn_mat)
        if ~cubo_conn_mat(i,j)
            continue
        end
        for k = 1:length(cubo_conn_mat)
            if ~cubo_conn_mat(j,k)
                continue
            end
            for l = 1:length(cubo_conn_mat)
                if ~cubo_conn_mat(k,l)
                    continue
                end
                for m = 1:length(cubo_conn_mat)
                    if ~cubo_conn_mat(l,m)
                        continue
                    end
                    for n = 1:length(cubo_conn_mat)
                        if ~cubo_conn_mat(m,n)
                            continue
                        end

                                
                                
                                if length([i,j,k,l,m,n]) ~= length(unique([i,j,k,l,m,n]))
                                    continue
                                end
                                
                                if cubo_conn_mat(i,j) & cubo_conn_mat(j,k) & cubo_conn_mat(k,l) & cubo_conn_mat(l,m) & cubo_conn_mat(m,n) & cubo_conn_mat(n,i)
                                    if isempty(cubo_vertex_sets_hex)
                                        cubo_vertex_sets_hex = [cubo_vertex_sets_hex;[i,j,k,l,m,n]];
                                    end
                                    if ~ismember(perms([i,j,k,l,m,n]),cubo_vertex_sets_hex,'rows')
                                        cubo_vertex_sets_hex = [cubo_vertex_sets_hex;[i,j,k,l,m,n]];
                                    end
                                end
                                
                                

                    end
                end
            end
        end
    end
end
% Octagons
for i = 1:length(cubo_conn_mat)
    for j = 1:length(cubo_conn_mat)
        if ~cubo_conn_mat(i,j)
            continue
        end
        for k = 1:length(cubo_conn_mat)
            if ~cubo_conn_mat(j,k)
                continue
            end
            for l = 1:length(cubo_conn_mat)
                if ~cubo_conn_mat(k,l)
                    continue
                end
                for m = 1:length(cubo_conn_mat)
                    if ~cubo_conn_mat(l,m)
                        continue
                    end
                    for n = 1:length(cubo_conn_mat)
                        if ~cubo_conn_mat(m,n)
                            continue
                        end
                        for o = 1:length(cubo_conn_mat)
                            if ~cubo_conn_mat(n,o)
                                continue
                            end
                            for p = 1:length(cubo_conn_mat)
                                if ~cubo_conn_mat(o,p)
                                    continue
                                end
                                
                                path = [i,j,k,l,m,n,o,p];
                                
                                % Check for any unwanted intermediate
                                % connections
                                intermediate_conn = 0;
                                for q = 1:length(path)
                                    shifted_path = circshift(path,2-q);
                                    for r = 4:length(path)
                                        if cubo_conn_mat(path(q),shifted_path(r))
                                            intermediate_conn = 1;
                                        end
                                    end
                                end
                                if intermediate_conn
                                    continue
                                end
                                
                                
                                if length(path) ~= length(unique(path))
                                    continue
                                end
                                

                                
                                if cubo_conn_mat(i,j) & cubo_conn_mat(j,k) & cubo_conn_mat(k,l) & cubo_conn_mat(l,m) & cubo_conn_mat(m,n) & cubo_conn_mat(n,o) & cubo_conn_mat(o,p) & cubo_conn_mat(p,i)
                                    if isempty(cubo_vertex_sets_oct)
                                        cubo_vertex_sets_oct = [cubo_vertex_sets_oct;path];
                                    end
                                    if ~ismember(perms(path),cubo_vertex_sets_oct,'rows')
                                        cubo_vertex_sets_oct = [cubo_vertex_sets_oct;path];
                                    end
                                end
                                
                                
                            end
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
% Squares
for i = 1:length(cubo_vertex_sets_square)
    fill3(...
        cubo(cubo_vertex_sets_square(i,:),1),...
        cubo(cubo_vertex_sets_square(i,:),2),...
        cubo(cubo_vertex_sets_square(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end
% Hexagons
for i = 1:length(cubo_vertex_sets_hex)
    fill3(...
        cubo(cubo_vertex_sets_hex(i,:),1),...
        cubo(cubo_vertex_sets_hex(i,:),2),...
        cubo(cubo_vertex_sets_hex(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end
% Octagons
for i = 1:size(cubo_vertex_sets_oct,1)
    fill3(...
        cubo(cubo_vertex_sets_oct(i,:),1),...
        cubo(cubo_vertex_sets_oct(i,:),2),...
        cubo(cubo_vertex_sets_oct(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)






