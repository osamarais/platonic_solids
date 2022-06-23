% Author: Osama M. Raisuddin

% This script creates 3D figures of the deltoidal icositetrahedron and its
% dual, the small rhombicuboctahedron



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



% Deltoidal icositetrahedron

% Write down the coordinates
sm_del_icos = [
     1, 0, 0;
    -1, 0, 0;
     
     0, 1, 0;
     0,-1, 0;
     
     0, 0, 1;
     0, 0,-1;
     
     
     0, 0.5*sqrt(2), 0.5*sqrt(2);
     0, 0.5*sqrt(2),-0.5*sqrt(2);
     0,-0.5*sqrt(2), 0.5*sqrt(2);
     0,-0.5*sqrt(2),-0.5*sqrt(2);
     
      0.5*sqrt(2), 0.5*sqrt(2),0;
      0.5*sqrt(2),-0.5*sqrt(2),0;
     -0.5*sqrt(2), 0.5*sqrt(2),0;
     -0.5*sqrt(2),-0.5*sqrt(2),0;
     
      0.5*sqrt(2), 0, 0.5*sqrt(2);
      0.5*sqrt(2), 0,-0.5*sqrt(2);
     -0.5*sqrt(2), 0, 0.5*sqrt(2);
     -0.5*sqrt(2), 0,-0.5*sqrt(2);
     
     
     (2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7;
     (2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7;
     (2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7;
     (2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7;
     
    -(2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7;
    -(2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7;
    -(2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7, (2*sqrt(2)+1)/7;
    -(2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7,-(2*sqrt(2)+1)/7;
    ];


% Generate the connectivity matrix based on the distances between the
% vertices
sm_del_icos_conn_mat = zeros(length(sm_del_icos),length(sm_del_icos));
for i = 1:length(sm_del_icos)
    for j = 1:length(sm_del_icos)
        if i==j
            continue
        end
        if (sum((sm_del_icos(i,:)-sm_del_icos(j,:)).^2,2) < 0.586) && (sum((sm_del_icos(i,:)-sm_del_icos(j,:)).^2,2) ~= 0)
            sm_del_icos_conn_mat(i,j) = 1;
        end
    end
end

% Generate the vertex sets using the connectivity matrix
sm_del_icos_vertex_sets = [];
for i = 1:length(sm_del_icos_conn_mat)
    for j = 1:length(sm_del_icos_conn_mat)
        for k = 1:length(sm_del_icos_conn_mat)
            for l = 1:length(sm_del_icos_conn_mat)
                
                if i==k || j==l
                    continue
                end
                
                if sm_del_icos_conn_mat(i,j) && sm_del_icos_conn_mat(j,k) && sm_del_icos_conn_mat(k,l) && sm_del_icos_conn_mat(l,i)
                    if isempty(sm_del_icos_vertex_sets)
                        sm_del_icos_vertex_sets = [sm_del_icos_vertex_sets;[i,j,k,l]];
                    end
                    if ~ismember(perms([i,j,k,l]),sm_del_icos_vertex_sets,'rows')
                        sm_del_icos_vertex_sets = [sm_del_icos_vertex_sets;[i,j,k,l]];
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
for i = 1:length(sm_del_icos_vertex_sets)
    fill3(...
        sm_del_icos(sm_del_icos_vertex_sets(i,:),1),...
        sm_del_icos(sm_del_icos_vertex_sets(i,:),2),...
        sm_del_icos(sm_del_icos_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)







% Small rhombicuboctahedron

% Write down the coordinates
sm_rhm_cub = [
     1, 1, (1+sqrt(2));
     1,-1, (1+sqrt(2));
    -1, 1, (1+sqrt(2));
    -1,-1, (1+sqrt(2));
    
     1, 1,-(1+sqrt(2));
     1,-1,-(1+sqrt(2));
    -1, 1,-(1+sqrt(2));
    -1,-1,-(1+sqrt(2));    
    
     1, (1+sqrt(2)), 1;
     1, (1+sqrt(2)),-1;
    -1, (1+sqrt(2)), 1;
    -1, (1+sqrt(2)),-1;
    
     1,-(1+sqrt(2)), 1;
     1,-(1+sqrt(2)),-1;
    -1,-(1+sqrt(2)), 1;
    -1,-(1+sqrt(2)),-1;
    
     (1+sqrt(2)), 1, 1;
     (1+sqrt(2)), 1,-1;
     (1+sqrt(2)),-1, 1;
     (1+sqrt(2)),-1,-1;
     
    -(1+sqrt(2)), 1, 1;
    -(1+sqrt(2)), 1,-1;
    -(1+sqrt(2)),-1, 1;
    -(1+sqrt(2)),-1,-1;
    ];

% Generate the connectivity matrix based on the distances between the
% vertices
sm_rhm_cub_conn_mat = zeros(length(sm_rhm_cub),length(sm_rhm_cub));
for i = 1:length(sm_rhm_cub)
    for j = 1:length(sm_rhm_cub)
        if i==j
            continue
        end
        if (sum((sm_rhm_cub(i,:)-sm_rhm_cub(j,:)).^2,2) <= 4) && (sum((sm_rhm_cub(i,:)-sm_rhm_cub(j,:)).^2,2) ~= 0)
            sm_rhm_cub_conn_mat(i,j) = 1;
        end
    end
end

% Generate the vertex sets using the connectivity matrix
sm_rhm_cub_vertex_sets_square = [];
sm_rhm_cub_vertex_sets_triangle = [];
% Squares
for i = 1:length(sm_rhm_cub_conn_mat)
    for j = 1:length(sm_rhm_cub_conn_mat)
        for k = 1:length(sm_rhm_cub_conn_mat)
            for l = 1:length(sm_rhm_cub_conn_mat)
                
                if i==k || j==l
                    continue
                end
                
                if sm_rhm_cub_conn_mat(i,j) && sm_rhm_cub_conn_mat(j,k) && sm_rhm_cub_conn_mat(k,l) && sm_rhm_cub_conn_mat(l,i)
                    if isempty(sm_rhm_cub_vertex_sets_square)
                        sm_rhm_cub_vertex_sets_square = [sm_rhm_cub_vertex_sets_square;[i,j,k,l]];
                    end
                    if ~ismember(perms([i,j,k,l]),sm_rhm_cub_vertex_sets_square,'rows')
                        sm_rhm_cub_vertex_sets_square = [sm_rhm_cub_vertex_sets_square;[i,j,k,l]];
                    end
                end
            end
        end
    end
end
% Triangles
for i = 1:length(sm_rhm_cub_conn_mat)
    for j = 1:length(sm_rhm_cub_conn_mat)
        for k = 1:length(sm_rhm_cub_conn_mat)
            if sm_rhm_cub_conn_mat(i,j) && sm_rhm_cub_conn_mat(j,k) && sm_rhm_cub_conn_mat(k,i)
                if isempty(sm_rhm_cub_vertex_sets_triangle)
                    sm_rhm_cub_vertex_sets_triangle = [sm_rhm_cub_vertex_sets_triangle;[i,j,k]];
                end
                if ~ismember(perms([i,j,k]),sm_rhm_cub_vertex_sets_triangle,'rows')
                    sm_rhm_cub_vertex_sets_triangle = [sm_rhm_cub_vertex_sets_triangle;[i,j,k]];
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
for i = 1:length(sm_rhm_cub_vertex_sets_square)
    fill3(...
        sm_rhm_cub(sm_rhm_cub_vertex_sets_square(i,:),1),...
        sm_rhm_cub(sm_rhm_cub_vertex_sets_square(i,:),2),...
        sm_rhm_cub(sm_rhm_cub_vertex_sets_square(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end
% Triangles
for i = 1:length(sm_rhm_cub_vertex_sets_triangle)
    fill3(...
        sm_rhm_cub(sm_rhm_cub_vertex_sets_triangle(i,:),1),...
        sm_rhm_cub(sm_rhm_cub_vertex_sets_triangle(i,:),2),...
        sm_rhm_cub(sm_rhm_cub_vertex_sets_triangle(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)

