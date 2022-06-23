% Author: Osama M. Raisuddin

% This script creates 3D figures of the 3N search vectors and its cosine
% measure vectors and compares the layout to the corresponding regular
% polyhedra

clc
clear
close all

color = [149, 121, 232]/255;
view_angle = [105,14];



%% 3N search and its cosine measure

% Write down the search directions
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







% Corresponding cosine measure vectors

% Generate the coordinates

% Write down the base, of whose permutations the vertices are formed
cm_base = [
     (sqrt(3) - sqrt(2)), (sqrt(2) - 1), 1;
     (sqrt(3) - sqrt(2)), (sqrt(2) - 1),-1;
     (sqrt(3) - sqrt(2)),-(sqrt(2) - 1), 1;
     (sqrt(3) - sqrt(2)),-(sqrt(2) - 1),-1;
    -(sqrt(3) - sqrt(2)), (sqrt(2) - 1), 1;
    -(sqrt(3) - sqrt(2)), (sqrt(2) - 1),-1;
    -(sqrt(3) - sqrt(2)),-(sqrt(2) - 1), 1;
    -(sqrt(3) - sqrt(2)),-(sqrt(2) - 1),-1;
     ];
% Make a single entry to avoid ismember() from failing
cm_vecs = [ (sqrt(3) - sqrt(2)), (sqrt(2) - 1), 1];
% Generate the remaining permutations
for i = 1:length(cm_base)
    myperms = perms(cm_base(i,:));
    for j = 1:length(myperms)
        if ~ismember(myperms(j,:),cm_vecs,'rows')
            cm_vecs = [cm_vecs; myperms(j,:)];
        end
    end
end








%% Deltoidal icositetrahedron and the great rhombicuboctahedron



% Deltoidal icositetrahedron (this forms the 'search directions')

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

% Normalize the vectors
for i = 1:length(sm_del_icos)
    sm_del_icos(i,:) = sm_del_icos(i,:)/norm(sm_del_icos(i,:));
end









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



%% The plots


% The geometric figures plot
figure
hold on
axis off
axis equal

% Scale down the size of the icositetrahedron
scale1 = 0.95;
% Scale up the size of the rhombicuboctahedron
scale2 = 2.5;
% Scale down the vectors of the rhombicuboctahedron
scale3 = scale2/0.9;
% Scale down the vectors of the rhombicuboctahedron
scale4 = scale2*0.7;
% Scale down the vectors of the icositetrahedron
scale5 = 0.6;

% Plot the edges of the deltoidal icositetrahedron
for i = 1:length(sm_del_icos_vertex_sets)
    plot3( scale1*sm_del_icos([sm_del_icos_vertex_sets(i,:) sm_del_icos_vertex_sets(i,1)],1) , scale1*sm_del_icos([sm_del_icos_vertex_sets(i,:) sm_del_icos_vertex_sets(i,1)],2) , scale1*sm_del_icos([sm_del_icos_vertex_sets(i,:) sm_del_icos_vertex_sets(i,1)],3), ':','Color','black','LineWidth',1)
end

% Plot the edges of the great rhombicuboctahedron
% Square
for i = 1:length(cubo_vertex_sets_square)
    plot3( scale2*cubo([cubo_vertex_sets_square(i,:) cubo_vertex_sets_square(i,1)],1) , scale2*cubo([cubo_vertex_sets_square(i,:) cubo_vertex_sets_square(i,1)],2) , scale2*cubo([cubo_vertex_sets_square(i,:) cubo_vertex_sets_square(i,1)],3), ':','Color',color,'LineWidth',1)
end
% Hexagon
for i = 1:length(cubo_vertex_sets_hex)
    plot3( scale2*cubo([cubo_vertex_sets_hex(i,:) cubo_vertex_sets_hex(i,1)],1) , scale2*cubo([cubo_vertex_sets_hex(i,:) cubo_vertex_sets_hex(i,1)],2) , scale2*cubo([cubo_vertex_sets_hex(i,:) cubo_vertex_sets_hex(i,1)],3), ':','Color',color,'LineWidth',1)
end
% Octagon
for i = 1:size(cubo_vertex_sets_oct,1)
    plot3( scale2*cubo([cubo_vertex_sets_oct(i,:) cubo_vertex_sets_oct(i,1)],1) , scale2*cubo([cubo_vertex_sets_oct(i,:) cubo_vertex_sets_oct(i,1)],2) , scale2*cubo([cubo_vertex_sets_oct(i,:) cubo_vertex_sets_oct(i,1)],3), ':','Color',color,'LineWidth',1)
end

% Plot vectors of the deltoidal icositetrahedron
quiver3(sm_del_icos(:,1)*scale5,sm_del_icos(:,2)*scale5,sm_del_icos(:,3)*scale5,sm_del_icos(:,1),sm_del_icos(:,2),sm_del_icos(:,3),'black','LineWidth',2)
% Plot vectors of the great rhombicuboctahedron
quiver3(cubo(:,1)*scale4,cubo(:,2)*scale4,cubo(:,3)*scale4,cubo(:,1)*scale3,cubo(:,2)*scale3,cubo(:,3)*scale3,'color',color,'LineWidth',2)

% Set viewing angle
view(view_angle)







% The 3N search plot
figure
hold on
axis off
axis equal


% Scale down the size of the 3N shape
scale1 = 0.95;
% Scale up the size of the CM shape
scale2 = 2.5;
% Scale down the vectors of the CM shape
scale3 = scale2/0.9;
% Scale down the vectors of the CM shape
scale4 = scale2*0.7;
% Scale down the vectors of the 3N shape
scale5 = 0.6;

% Plot the edges of the 3N shape
for i = 1:length(sm_del_icos_vertex_sets)
    plot3( scale1*search_vecs([sm_del_icos_vertex_sets(i,:) sm_del_icos_vertex_sets(i,1)],1) , scale1*search_vecs([sm_del_icos_vertex_sets(i,:) sm_del_icos_vertex_sets(i,1)],2) , scale1*search_vecs([sm_del_icos_vertex_sets(i,:) sm_del_icos_vertex_sets(i,1)],3), ':','Color','black','LineWidth',1)
end

% Plot the edges of the CM shape
% Square
for i = 1:length(cubo_vertex_sets_square)
    plot3( scale2*cm_vecs([cubo_vertex_sets_square(i,:) cubo_vertex_sets_square(i,1)],1) , scale2*cm_vecs([cubo_vertex_sets_square(i,:) cubo_vertex_sets_square(i,1)],2) , scale2*cm_vecs([cubo_vertex_sets_square(i,:) cubo_vertex_sets_square(i,1)],3), ':','Color',color,'LineWidth',1)
end
% Hexagon
for i = 1:length(cubo_vertex_sets_hex)
    plot3( scale2*cm_vecs([cubo_vertex_sets_hex(i,:) cubo_vertex_sets_hex(i,1)],1) , scale2*cm_vecs([cubo_vertex_sets_hex(i,:) cubo_vertex_sets_hex(i,1)],2) , scale2*cm_vecs([cubo_vertex_sets_hex(i,:) cubo_vertex_sets_hex(i,1)],3), ':','Color',color,'LineWidth',1)
end
% Octagon
for i = 1:size(cubo_vertex_sets_oct,1)
    plot3( scale2*cm_vecs([cubo_vertex_sets_oct(i,:) cubo_vertex_sets_oct(i,1)],1) , scale2*cm_vecs([cubo_vertex_sets_oct(i,:) cubo_vertex_sets_oct(i,1)],2) , scale2*cm_vecs([cubo_vertex_sets_oct(i,:) cubo_vertex_sets_oct(i,1)],3), ':','Color',color,'LineWidth',1)
end

% Plot vectors of the search vectors
quiver3(search_vecs(:,1)*scale5,search_vecs(:,2)*scale5,search_vecs(:,3)*scale5,search_vecs(:,1),search_vecs(:,2),search_vecs(:,3),'black','LineWidth',2)
% Plot vectors of the cosine measure vectors
quiver3(cm_vecs(:,1)*scale4,cm_vecs(:,2)*scale4,cm_vecs(:,3)*scale4,cm_vecs(:,1)*scale3,cm_vecs(:,2)*scale3,cm_vecs(:,3)*scale3,'color',color,'LineWidth',2)

% Set viewing angle
view(view_angle)

