% Author: Osama M. Raisuddin

% This script creates 3D figures of the searches done with regular convex
% polyhedra and their cosine measures and duality property

close all
clear
clc


% Use the high quality renderer to ensure SVG plots are created
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')





% Choose the
% color,
% and viewing angle
% for the figures
color = [149, 121, 232]/255;
view_angle = [120,30];




%% Nelder Mead

% Generate the search directions for Nelder Mead search and the vertex sets
% for the corresponding polyhedron, the tetrahedron
search_dir = [1, 1, 1;1, -1, -1;-1, 1, -1 ;-1, -1,  1]/sqrt(3);
search_dir_vertex_sets = [[1,2,3];[1,2,4];[1,3,4];[2,3,4]];

% Generate the cosine measure directions for Nelder Mead search and the 
% vertex sets for the corresponding polyhedron, the tetrahedron
cm_dir = [-1, -1, -1; -1, 1, 1; 1, -1, 1; 1,  1, -1]/sqrt(3);
cm_dir_vertex_sets = [[1,2,3];[1,2,4];[1,3,4];[2,3,4]];

% Plot the polyhera and vectors
figure
hold on
axis off
axis equal

quiver3(zeros(length(search_dir),1),zeros(length(search_dir),1),zeros(length(search_dir),1),search_dir(:,1),search_dir(:,2),search_dir(:,3),'black','LineWidth',2)
quiver3(zeros(length(cm_dir),1),zeros(length(cm_dir),1),zeros(length(cm_dir),1),cm_dir(:,1),cm_dir(:,2),cm_dir(:,3),'color',color,'LineWidth',2)

for i = 1:length(search_dir_vertex_sets)
    plot3( search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],1) , search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],2) , search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],3), ':','Color','black')
    
end
for i = 1:length(cm_dir_vertex_sets)
    plot3( cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],1) , cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],2) , cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],3), ':','Color',color)
    
end
view(view_angle)




%%  Pattern Search
search_dir = [1,0,0;-1,0,0; 0,1,0;0,-1,0; 0,0,1;0,0,-1];
search_dir_vertex_sets = [1,3,5; 1,3,6; 1,4,5; 1,4,6; 2,3,5; 2,3,6; 2,4,5; 2,4,6;];


cm_dir = [
     1, 1, 1;
     1, 1,-1;
     1,-1, 1;
     1,-1,-1;
    -1, 1, 1;
    -1, 1,-1;
    -1,-1, 1;
    -1,-1,-1;
     ]/sqrt(3);
cm_dir_vertex_sets = [
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

 
 
quiver3(zeros(length(search_dir),1),zeros(length(search_dir),1),zeros(length(search_dir),1),search_dir(:,1),search_dir(:,2),search_dir(:,3),'black','LineWidth',2)
quiver3(zeros(length(cm_dir),1),zeros(length(cm_dir),1),zeros(length(cm_dir),1),cm_dir(:,1),cm_dir(:,2),cm_dir(:,3),'color',color,'LineWidth',2)
for i = 1:length(search_dir_vertex_sets)
    plot3( search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],1) , search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],2) , search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],3), ':','Color','black')
    
end
for i = 1:length(cm_dir_vertex_sets)
    plot3( cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],1) , cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],2) , cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],3), ':','Color',color)
    
end
view(view_angle)




%%  Hypercube Search
search_dir = [
     1, 1, 1;
     1, 1,-1;
     1,-1, 1;
     1,-1,-1;
    -1, 1, 1;
    -1, 1,-1;
    -1,-1, 1;
    -1,-1,-1;
     ]/sqrt(3);
search_dir_vertex_sets = [
    1,2,4,3;
    1,5,7,3;
    1,2,6,5;
    8,7,5,6;
    8,4,3,7;
    8,4,2,6
    ]; 
 
 
 
cm_dir = [1,0,0;-1,0,0; 0,1,0;0,-1,0; 0,0,1;0,0,-1];
cm_dir_vertex_sets = [1,3,5; 1,3,6; 1,4,5; 1,4,6; 2,3,5; 2,3,6; 2,4,5; 2,4,6;];
 

figure
hold on
axis off
axis equal

 
 
quiver3(zeros(length(search_dir),1),zeros(length(search_dir),1),zeros(length(search_dir),1),search_dir(:,1),search_dir(:,2),search_dir(:,3),'black','LineWidth',2)
quiver3(zeros(length(cm_dir),1),zeros(length(cm_dir),1),zeros(length(cm_dir),1),cm_dir(:,1),cm_dir(:,2),cm_dir(:,3),'color',color,'LineWidth',2)
for i = 1:length(search_dir_vertex_sets)
    plot3( search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],1) , search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],2) , search_dir([search_dir_vertex_sets(i,:) search_dir_vertex_sets(i,1)],3), ':','Color','black')
    
end
for i = 1:length(cm_dir_vertex_sets)
    plot3( cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],1) , cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],2) , cm_dir([cm_dir_vertex_sets(i,:) cm_dir_vertex_sets(i,1)],3), ':','Color',color)
    
end

view(view_angle)


