clc
clear
close all

% Use high quality renderer if SVG plots are to be created
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')



% The height of the plane showing the movement of the algorithms
h = 15;



%% Pattern Search


% Resolution of the surface of the functional. Going above 800 stalls the
% code
surface_resolution = 100;


% Evaluate the peaks function and get the XYZ values
[X,Y,Z] = peaks(surface_resolution);

% Do this to crop the surface in a circular shape
Z(X.^2+Y.^2>9) = NaN;


% Generate the surface
figure
hold on
% Custom color map
% map = [58, 32, 140
% 202, 189, 242]/255;
%
% map_values = 10000;
% map = map(1,:) + kron(map(2,:) - map(1,:),(0:map_values)'/map_values);
%
% colormap(map)
surf(X,Y,Z)
shading interp
axis off



% Initial guess and search size
u = [0.7,0.4];
delta = 1;


% Perform a few steps of pattern search and show it on the surface
for i = 1:6
    
    [x,y] = meshgrid(u(1)+delta*(-1:1),u(2)+delta*(-1:1));
    
    z = peaks(x,y);
    
    % The lines of the patter
    plot3(x(2,[1,3]),y(2,[1,3]),h*ones(1,2),'black','LineWidth',2);
    plot3(x([1,3],2),y([1,3],2),h*ones(1,2),'black','LineWidth',2);
    
    % The lines from the plane to the surface
    plot3([x(2,2) x(2,2)],[y(2,2) y(2,2)],[h peaks(x(2,2),y(2,2))+1],'black','LineWidth',1,'LineStyle',':');
    
    % The dots of the pattern
    scatter3([x(2,[1,3]),x([1,3],2)'],[y(2,[1,3]),y([1,3],2)'],ones(1,4)*h,ones(1,4)*50,'o','filled','red')
    % The dots on the surface (best guesses)
    scatter3(x(2,2),y(2,2),peaks(x(2,2),y(2,2))+1,100,'o','filled','black')

    % The best best guesses in the plane
    scatter3(x(2,2),y(2,2),h,120,'o','black','LineWidth', 1.5)
    
    
    
    % Flag for failed iteration
    flag = false;
    
    for j = 1:3
        if z(2,2)>z(2,j)
            u(1) = x(2,j);
            u(2) = y(2,j);
            z(2,2) = z(2,j);
            flag = true;
        end
        if z(2,2)>z(j,2)
            u(1) = x(j,2);
            u(2) = y(j,2);
            z(2,2) = z(j,2);
            flag = true;
        end
    end
    
    if ~flag
        delta = delta/2;
    end
    
    
    
end


% Set the z limits to make the plot look good
zlim([-8 h+5]);
% Set the viewing angle
view([-60 69])

% Save the plots as an image file
%saveas(gcf,'pattern_search.svg')




%% Nelder Mead

% Run the Nelder Mead code and paste the output here
nelder_mead_results = [0.7000    0.4000
    1.7000    0.4000
    0.7000    1.4000
    
    0.7000    0.4000
    1.7000   -0.6000
    1.7000    0.4000
    
    0.2000   -1.1000
    0.7000    0.4000
    1.7000   -0.6000
    
    0.2000   -1.1000
    -0.8000   -0.1000
    0.7000    0.4000
    
    0.2000   -1.1000
    0.4500   -0.3500
    -0.3000   -0.6000
    
    0.2000   -1.1000
    0.6375   -0.7875
    0.4500   -0.3500
    
    0.3875   -1.5375
    0.2000   -1.1000
    0.6375   -0.7875
    
    0.3875   -1.5375
    -0.0500   -1.8500
    0.2000   -1.1000];




% Generate the surface
figure
hold on
% Custom color map
% map = [58, 32, 140
% 202, 189, 242]/255;
%
% map_values = 10000;
% map = map(1,:) + kron(map(2,:) - map(1,:),(0:map_values)'/map_values);
%
% colormap(map)
surf(X,Y,Z)
shading interp
axis off



% Plot the progress of the algorithm on the surface
for i = 1:3:length(nelder_mead_results)
    
    points = nelder_mead_results(i:(i+2),:);
    
    z = peaks(points(:,1),points(:,2));
    
    [~,index] = min(z);
    
    
    % Triangle vertices
    scatter3(points(:,1),points(:,2),ones(3,1)*h, ones(3,1)*50,'o','filled','red')
    plot3([points(:,1);points(1,1)],[points(:,2);points(1,2)],ones(4,1)*h,'black','LineWidth',2);
    
    if i == 1
        % The ball on the surface
        scatter3(points(index,1),points(index,2),z(index)+0.5,100,'o','filled','black')
        % Best points of the triangles
        scatter3(points(index,1),points(index,2),h,120,'o','black','LineWidth', 1.5)
        % Dotted line b/w top and surface
        plot3([points(index,1) points(index,1)],[points(index,2) points(index,2)],[h z(index)+1],'black','LineWidth',1,'LineStyle',':');
    
    else
        for j = 1:3
            if ~ ismember(points(j,:),nelder_mead_results((i-3):(i-1),:),'rows')
                scatter3(points(j,1),points(j,2),h,120,'o','black','LineWidth', 1.5)
                scatter3(points(j,1),points(j,2),z(j)+0.5,100,'o','filled','black')
                plot3([points(j,1) points(j,1)],[points(j,2) points(j,2)],[h z(j)+0.5],'black','LineWidth',1,'LineStyle',':');
            end
        end
    end
    
    
    
    
end



% Set the z limits to make the plot look good
zlim([-8 h+5]);
% Set the viewing angle
view([-60 69])

% Save the plots as an image file
%saveas(gcf,'nelder_mead.svg')




% Ax = gca;
% Ax.ZAxis.Visible = 'off';
% Ax.ZGrid = 'off';
% Ax.Color = 'none';