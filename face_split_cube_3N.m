% Author: Osama M. Raisuddin

% This script creates a 3D figure of the cube with its faces split for the
% 3N search method


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





% Cube with split faces

% Write down the coordinates
search_vecs = [1,0,0;
    -1,0,0;
    0,1,0;
    0,-1,0;
    0,0,1;
    0,0,-1;
    0,1,1;
    0,1,-1;
    0,-1,1;
    0,-1,-1;
    1,1,0;
    1,-1,0;
    -1,1,0;
    -1,-1,0;
    1,0,1;
    1,0,-1;
    -1,0,1;
    -1,0,-1;
    1,1,1;
    1,1,-1;
    1,-1,1;
    1,-1,-1;
    -1,1,1;
    -1,1,-1;
    -1,-1,1;
    -1,-1,-1];

figure
hold on
axis off
axis equal

% The vertex sets are the same as the small deltoidal icositetrahedron!
% The code is therefore borrowed

% Generate the connectivity matrix based on the distances between the
% vertices
sm_del_icos_conn_mat = zeros(length(search_vecs),length(search_vecs));
for i = 1:length(search_vecs)
    for j = 1:length(search_vecs)
        if i==j
            continue
        end
        if (sum((search_vecs(i,:)-search_vecs(j,:)).^2,2) < 1.1) && (sum((search_vecs(i,:)-search_vecs(j,:)).^2,2) ~= 0)
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

% Generate faces of the solids
for i = 1:length(sm_del_icos_vertex_sets)
    fill3(...
        search_vecs(sm_del_icos_vertex_sets(i,:),1),...
        search_vecs(sm_del_icos_vertex_sets(i,:),2),...
        search_vecs(sm_del_icos_vertex_sets(i,:),3),...
        color,...
        'FaceAlpha',alpha)
end

% Set the viewing angle
view(view_angle)

% The figure does not show up properly as an SVG
% Export it as a high quality PNG instead if necessary
% exportgraphics(gcf,'face_split_cube.png','Resolution',1000)
