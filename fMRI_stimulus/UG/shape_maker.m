% making of shapes for ultimatum game
% Author: Michael Giffin, January, 2018

%% Create  Filled Circle
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);
% 
% xlim([1 300]);
% ylim([1 300]);
h = patch;
axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,120);
h.XData = 256 + 128*cos(t);
h.YData = 256 + 128*sin(t);
h.FaceVertexCData = h.XData;
% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/circle_filled.png')

%% Create Empty Circle
figure('Position',[100 100 512 512],'Color','white')
axes('Position',[0 0 1 1])
xlim([1 512]);
ylim([1 512]);

h = line;
axis off
colormap(gray)

% Set Coordinates
t = linspace(0,2*pi,120);
h.XData = 256 + 128*cos(t);
h.YData = 256 + 128*sin(t);
set(h(1),'linewidth',10);

% Save Image
% img = getframe(gcf);
% imwrite(img.cdata,'images/circle_frame.png')

%% Create filled Square
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [100 300 300 100 100];
y = [100 100 300 300 100];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/square_filled.png')


%% Create empty Square
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

% h = patch;
axis off
colormap(gray)

x = [100 300 300 100 100];
y = [100 100 300 300 100];

patch(x, y, 'white', 'linewidth', 10);

% h.XData = x;
% h.YData = y;
% h.FaceVertexCData = h.XData;

% 
% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/square_frame.png')


%% Create filled Triangle
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [100 300 200 100];
y = [100 100 300 100];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/triangle_filled.png')


%% Create empty Triangle
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [100 300 200 100];
y = [100 100 300 100];

patch(x, y, 'white', 'linewidth', 10);

% h.XData = x;
% h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/triangle_frame.png')


%% Create filled Diamond
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

hyp = 200/(sqrt(2)); % hypotenuse or diagonal of square to keep size the same
% hyp = 100;

x = [55 55+hyp 55+hyp+hyp 55+hyp 55];
y = [55+hyp 55 55+hyp 55+hyp+hyp 55+hyp];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/diamond_filled.png')


%% Create empty Diamond
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

hyp = 200/(sqrt(2)); % hypotenuse or diagonal of square to keep size the same

x = [55 55+hyp 55+hyp+hyp 55+hyp 55];
y = [55+hyp 55 55+hyp 55+hyp+hyp 55+hyp];

patch(x, y, 'white', 'linewidth', 10);


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/diamond_frame.png')


%% Create filled trapazoid
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [160 240 330 70 160];
y = [100 100 300 300 100];

h.XData = x;
h.YData = y;
% h.FaceVertexCData = h.XData;

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/trapazoid_filled.png')


%% Create empty trapazoid
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

h = patch;
axis off
colormap(gray)

x = [160 240 330 70 160];
y = [100 100 300 300 100];

patch(x, y, 'white', 'linewidth', 10);

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/trapazoid_frame.png')

%% Create filled star
% I had to fill it with a circle, it took longer than it should have to get
% this right
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);

xCenter = 200;
yCenter = 200;

h = patch;
hold on
% to fill interior with circle
% Set Coordinates
t = linspace(0,2*pi,120);
h.XData = xCenter + 57*cos(t);
h.YData = yCenter + 57*sin(t);
h.FaceVertexCData = h.XData;


% Define parameters
numberOfPoints = 5;
rotationAngle = 0;
xCenter = 200;
yCenter = 200;
% Determine the angles that the arm tips are at
theta = (0 : (numberOfPoints-1)/numberOfPoints*pi : (numberOfPoints-1)*pi) + rotationAngle;
% Define distance from the arm tip to the center of star.
amplitude = 150;
% Get x and y coordinates of the arm tips.
x = amplitude .* cos(theta) + xCenter;
y = amplitude .* sin(theta) + yCenter;


h = patch;
axis off
colormap(gray)


h.XData = x;
h.YData = y;
h.FaceVertexCData = h.XData;

% % the distance between two points on the star
% ptdis = sqrt(((y(4) - y(1))^2) + ((x(4) - x(1))^2));
% 
% % the angles of the resulting triangle are 36, 36, and 108
% % the length of the edge of one of the stars triangles is
% triedge = (ptdis * sind(36))/(sind(108));
% 
% % and the length of the interior of the star's pentagon
% pentedge = (x(1) - x(2)) - (2*triedge);


% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/star_filled.png')


%% Create empty star
figure('Position',[0 0 550 550],'Color','white')
axes('Position',[0 0 1 1])
xlim([0 400]);
ylim([0 400]);


% Define parameters
numberOfPoints = 5;
rotationAngle = 0;
xCenter = 200;
yCenter = 200;
% Determine the angles that the arm tips are at
theta = (0 : (numberOfPoints-1)/numberOfPoints*pi : (numberOfPoints-1)*pi) + rotationAngle;
% Define distance from the arm tip to the center of star.
amplitude = 150;
% Get x and y coordinates of the arm tips.
x = amplitude .* cos(theta) + xCenter;
y = amplitude .* sin(theta) + yCenter;

h = patch;
axis off
colormap(gray)

patch(x, y, 'white', 'linewidth', 10);

% Save Image
img = getframe(gcf);
imwrite(img.cdata,'images/star_frame.png')

