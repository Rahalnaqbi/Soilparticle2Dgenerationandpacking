%Generating 2D grain shapes based on stochastic analysis of diff samples

%1: 2D projections and required values of contour 
%Extract aspect ratio and normalized Fourier coefficients for all contours
%/////////////////////////////////////////////////////////////////////////

%Initalize variables 
As = [];
folderPath = uigetdir; 
stlFiles = dir(fullfile(folderPath, '*.stl'));

for k = 1:length(stlFiles)
    fullFilePath = fullfile(folderPath, stlFiles(k).name);
end

excelFile = dir(fullfile(folderPath, '*.xlsx'));
exfilename = excelFile.name;
stlFileColumnName = 'Var47'; 
ASColumnName = 'Var43';

excelData = readtable(exfilename, 'Sheet', '10084', 'VariableNamingRule', 'preserve');
aan = zeros(20,4700);
bbn = zeros(20,4700);
aam = zeros(100,4700);
bbm = zeros(100,4700);


for i = 1:4700
stlFile = fullfile(folderPath, stlFiles(i).name);
%ensure that coordinates of 2d contour are found 
model = stlread(stlFile);
vertices = model.Points;
X = vertices(:,1);
Y = vertices(:,2);
Z = vertices(:,3);

desiredZ = median(Z); % Replace with your target Z height

% Find the indices of points around the desired Z level
tolerance = 0.1; % Set a tolerance for matching Z values
indices = abs(Z - desiredZ) < tolerance;

% Extract corresponding X and Y for that Z value
x = X(indices);
y = Y(indices);

% Plot the outer contour
points = [x, y];
%plot(x, y, 'o');

% Initialize an array to hold connections
distances = pdist2(points, points);
distances(logical(eye(size(distances)))) = Inf;

[~, fileName, fileExt] = fileparts(stlFile);
fullFileName = [fileName, fileExt];
disp(fullFileName);
rowIdx = find(strcmp(excelData.(stlFileColumnName), fullFileName));

STL_AS = excelData.(ASColumnName)(rowIdx);
As(i)= STL_AS;

centroidX = mean(x);
centroidY = mean(y);

% Shift the contour points to center them at the origin
xShifted = x - centroidX;
yShifted = y - centroidY;

% Convert the contour points to polar coordinates (r, theta)
[thetaa, r] = cart2pol(xShifted, yShifted);
[thetaSorted, sortIdx] = sort(thetaa);
rSorted = r(sortIdx);
N = length(rSorted);  % Number of points
fourierCoeffs = fft(rSorted);
max_coeff1 = 20;
max_coeff2 = 100; 
a_n = real(fourierCoeffs(1:min(max_coeff1, length(fourierCoeffs))));
b_n = imag(fourierCoeffs(1:min(max_coeff1, length(fourierCoeffs))));

if length(a_n) < max_coeff1
    a_n = [a_n; zeros(max_coeff1 - length(a_n), 1)];
end

if length(b_n) < max_coeff1
    b_n = [b_n; zeros(max_coeff1 - length(b_n), 1)];
end

if length(fourierCoeffs) >= 20
a_m = real(fourierCoeffs(max_coeff1+1:min(max_coeff1+max_coeff2, length(fourierCoeffs))));
b_m = imag(fourierCoeffs(max_coeff1+1:min(max_coeff1+max_coeff2, length(fourierCoeffs))));
if length(a_m) < max_coeff2
    a_m = [a_m; zeros(max_coeff2 - length(a_m), 1)];
end

if length(b_m) < max_coeff2
    b_m = [b_m; zeros(max_coeff2 - length(b_m), 1)];
end
else
    a_m = zeros(100,1);
    b_m = zeros(100,1);
end

a_n = a_n / max(a_n);
b_n = b_n / max(b_n);
a_m = a_m / max(a_m);
b_m = b_m / max(b_m);
aan(:,i) = a_n;
bbn(:,i) = b_n;
aam(:,i) = a_m;
bbm(:,i) = b_m;

end

%/////////////////////////////////////////////////////////////////////////
aam = aam(:, ~any(isnan(aam), 1));
bbm = bbm(:, ~any(isnan(bbm), 1));
aan = aan(:, ~any(isnan(aan), 1));
bbn = bbn(:, ~any(isnan(bbn), 1));
aam = abs(aam);
aan = abs(aan);
bbm = abs(bbm);
bbn = abs(bbn);

%2: Average required function values of all particle 2D contours
%/////////////////////////////////////////////////////////////////////////
avg_an = mean(aan, 2);
avg_bn = mean(bbn, 2);
avg_am = mean(aam, 2);
avg_bm = mean(bbm, 2);
avg_As = mean(As);

%/////////////////////////////////////////////////////////////////////////

%3: Constructing Fourier spectrum of soil sample
%/////////////////////////////////////////////////////////////////////////

% Define parameters for the base ellipse (aspect ratio)
shortAxis = 1;  
aspectRatio = avg_As;

% Generate the base ellipse with aspect ratio
numPoints = 1000;
theta = linspace(0, 2*pi, numPoints);  
AS = @(theta) (aspectRatio * shortAxis) ./ sqrt((shortAxis * cos(theta)).^2 + (aspectRatio * sin(theta)).^2);

r = @(theta) AS(theta);

% Add angularity to the base ellipse using average Fourier coefficients (20-order)
numAngularityTerms = length(avg_an);
for n = 1:numAngularityTerms
    r = @(theta) r(theta) + avg_an(n) * cos(n * theta) + avg_bn(n) * sin(n * theta);  % Angularity terms
end

% Add surface texture using average Fourier coefficients (100-order)
numTextureTerms = length(avg_am);
for m = 1:numTextureTerms
    r = @(theta) r(theta) + avg_am(m) * cos((m + numAngularityTerms) * theta) + avg_bm(m) * sin((m + numAngularityTerms) * theta);  % Surface texture terms
end

%Constructing normalized amplitude spectrum 
A_n = sqrt(avg_an.^2 + avg_bn.^2);

A_m = sqrt(avg_am.^2 + avg_bm.^2);

normalizedAmplitudeSpectruman = A_n / aspectRatio;
normalizedAmplitudeSpectrumam = A_m / aspectRatio;

%figure;
%stem(1:length(normalizedAmplitudeSpectruman), normalizedAmplitudeSpectruman, 'LineWidth', 2);
%xlabel('Mode Number');
%ylabel('Normalized Amplitude (Angularity)');
%title('Averaged Fourier Descriptors Values');
%grid on;

%figure;
%stem(1:length(normalizedAmplitudeSpectrumam), normalizedAmplitudeSpectrumam, 'LineWidth', 2);
%xlabel('Mode Number');
%ylabel('Normalized Amplitude (Surface texture)');
%title('Averaged Fourier Descriptors Values');
%grid on;

%/////////////////////////////////////////////////////////////////////////

%4: Defining packing and generation parameters  
%/////////////////////////////////////////////////////////////////////////

%Define number of generated particles needed
Numparticles = 250;

%Generate particles required 
r = zeros(numPoints,Numparticles);
theta = linspace(0, 2*pi, numPoints);

for i = 1:Numparticles

delta_n = -pi + 2*pi*rand(20, 1);
A_nn = normalizedAmplitudeSpectruman .* cos(delta_n);
B_nn = normalizedAmplitudeSpectruman .* sin(delta_n);

delta_m = -pi + 2*pi*rand(100, 1);
A_mm = normalizedAmplitudeSpectrumam .* cos(delta_m);
B_mm = normalizedAmplitudeSpectrumam .* sin(delta_m);

firstv = [];
secondv = [];
thirdv = [];
%Generating particle points
for t = 1:numPoints
    firstv(t) = aspectRatio;
for a = 1:20
    secondv(t) = A_nn(a) * cos((a) * theta(t)) + B_nn(a) * sin((a) * theta(t));
end
for b = 1:100
    thirdv(t) = A_mm(b) * cos((b) * theta(t)) + B_mm(b) * sin((b) * theta(t));
end
r(t,i) = firstv(t) + secondv(t) + thirdv(t);
end


end

r_values = r(:,1);
r_values2 = r(:,100);
figure;
polarplot(theta,r_values);
figure;
polarplot(theta,r_values2);
return;
%Packing generated particles 
containerWidth = 100;
containerHeight = 100;
targetSolidFraction = 0.9; 
minDist = 5;
particleCenters = zeros(Numparticles, 2);
% Initial random points for particle centers 
for i = 1:Numparticles
    validPosition = false;
    while ~validPosition
        newCenter = [containerWidth * rand, containerHeight * rand];
        distances = sqrt(sum((particleCenters(1:i-1, :) - newCenter).^2, 2));
        if all(distances >= minDist)
            particleCenters(i, :) = newCenter;
            validPosition = true;
        end
    end
end

surfaceareacolumn = 'Var2';
surfaceareas = excelData.(surfaceareacolumn);
meanSurfaceArea = mean(surfaceareas);  
stdSurfaceArea = std(surfaceareas); 

%Target size distribution 
targetSizes = normrnd(meanSurfaceArea, stdSurfaceArea, [Numparticles, 1]);
normtargetSizes = targetSizes / max(targetSizes);
%Voronoi Cells

%Initial cells areas
voronoiAreas = computeVoronoiAreas(particleCenters, containerWidth, containerHeight);

%Fitting cells to meet target size distribution using IMC algorithm 
maxIterations = 5000;
errorThreshold = 0.1;
currentError = inf;

for iter = 1:maxIterations
    [voronoiAreas, triangleX, triangleY] = computeVoronoiAreas(particleCenters, containerWidth, containerHeight);
    differenceMatrix = zeros(Numparticles, Numparticles);
    normvoronoiAreas = voronoiAreas / max(voronoiAreas);
for i = 1:Numparticles
    for j = 1:Numparticles
        differenceMatrix(i, j) = (normvoronoiAreas(i) - normtargetSizes(j))^2;
    end
end
[assignments, currentError] = matchpairs(differenceMatrix, 1e10);

    if currentError < errorThreshold
        disp('Converged');
        break;
    end
 
    randIndex = randi(Numparticles);
    perturbation = 0.1 * randn(1, 2);  
    newCenters = particleCenters;
    newCenters(randIndex, :) = particleCenters(randIndex, :) + perturbation;
    
    [newVoronoiAreas,~,~] = computeVoronoiAreas(newCenters, containerWidth, containerHeight);
    normnewVoronoiAreas = newVoronoiAreas / max(newVoronoiAreas);
    for i = 1:Numparticles
    for j = 1:Numparticles
        differenceMatrix2(i, j) = (normnewVoronoiAreas(i) - normtargetSizes(j))^2;
    end
end
[assignments, newError] = matchpairs(differenceMatrix2, 1e10);
    
    if newError < currentError
        particleCenters = newCenters;
        currentError = newError;
    end
end

figure;
[vx, vy] = voronoi(particleCenters(:,1), particleCenters(:,2));
numCenters = size(particleCenters, 1);
voronoiCells = cell(numCenters, 1);
uniqueCells = false(numCenters, 1);
% Delaunay triangulation for indexing
dt = delaunayTriangulation(particleCenters);

% Voronoi cell indices for each particle center
for i = 1:numCenters
    
    voronoiIndices = dt.ConnectivityList(dt.ConnectivityList(:, 1) == i, :);
    if ~uniqueCells(i)
        
        voronoiCells{i} = particleCenters(i, :);
        uniqueCells(i) = true; 
    end
    % Collecting particle centers in the corresponding Voronoi cell
    for j = 1:size(voronoiIndices, 1)
        idx = voronoiIndices(j, :);
        
        % Only keeping the first particle center for each Voronoi cell
        for k = 1:length(idx)
            if ~uniqueCells(idx(k)) 
                voronoiCells{idx(k)} = particleCenters(idx(k), :); 
                uniqueCells(idx(k)) = true; 
            end
        end
    end
end

areasfinal = voronoiAreas(uniqueCells);
cellX = triangleX(:,uniqueCells);
cellY = triangleY(:,uniqueCells);

plot(vx, vy, '-');
hold on;
for i = 1:length(voronoiCells)
    if ~isempty(voronoiCells{i}) 
        plot(voronoiCells{i}(1), voronoiCells{i}(2), 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
end
xlim([0 100]); 
ylim([0 100]);
title('Final Voronoi Tessellation Optimized for Target Size Distribution');

theta = theta(:);
%Fitting generated particles into cells
optimizedCenters = voronoiCells; 
optimizedCenterX = [];
optimizedCenterY = [];
particleX = zeros(numPoints,length(voronoiCells));
particleY = zeros(numPoints,length(voronoiCells));
for i = 1:length(voronoiCells)
    targetArea = areasfinal(i);
    particleX(:,i) = r(:, i) .* cos(theta) + voronoiCells{i}(1);
    particleY (:,i) =  r(:, i) .* sin(theta) + voronoiCells{i}(2);
    particleX2 = voronoiCells{i}(1);
    particleY2 = voronoiCells{i}(2);
    cellX2 = cellX (:,i);
    cellY2 = cellY(:,i);
    initialParams = [0, 0, 1];
    r2 = r(:, i);
    options = optimoptions('fmincon', 'Display', 'off');
    optimizedParams = fmincon(@(params) objective(params, r2, theta, cellX2, cellY2, particleX2, particleY2, targetArea, targetSolidFraction), ...
                              initialParams, [], [], [], [], [-inf, -inf, 0.1], [inf, inf, 2], [], options);

    dx = optimizedParams(1);
    dy = optimizedParams(2);
    scaleFactor = optimizedParams(3);

    optimizedCenterX(i) = optimizedCenters{i}(1) + dx;
    optimizedCenterY(i) = optimizedCenters{i}(2) + dy;
    r(:, i) = scaleFactor * r(:, i);

    particleX (:,i) = r(:, i) .* cos(theta) + optimizedCenterX(i);
    particleY (:,i) = r(:, i) .* sin(theta) + optimizedCenterY(i);

end

%Final packing 
figure;
hold on;
set(gca, 'color', 'k');  
for i = 1:length(voronoiCells)
    
    particleX (:,i) = r(:, i) .* cos(theta) + optimizedCenterX(i);
    particleY (:,i) = r(:, i) .* sin(theta) + optimizedCenterY(i);
    
   
    fill(particleX, particleY, 'w');  
end
axis equal;
title('Final Packed Particles (without Voronoi cells)');
hold off;

height = 0.5; 
numParticles = length(voronoiCells); 

allVertices = [];
faces = [];

function [areas, triangleX, triangleY] = computeVoronoiAreas(centers, width, height)
    [vx, vy] = voronoi(centers(:,1), centers(:,2));
    validIndices = (vx >= 0 & vx <= width) & (vy >= 0 & vy <= height);
    vx = vx(validIndices);
    vy = vy(validIndices);
    [polygons, ~] = voronoiDiagram(delaunayTriangulation(centers));
    finitepolygons = all(isfinite(polygons), 2);
    filteredpolygons = polygons(finitepolygons, :);
    filteredpolygons = abs(filteredpolygons);
    isWithinBounds = (filteredpolygons(:, 1) >= 0 & filteredpolygons(:, 1) <= width) & ...
                 (filteredpolygons(:, 2) >= 0 & filteredpolygons(:, 2) <= height);
    filteredPolygons = filteredpolygons(isWithinBounds, :);
    numCells = length(filteredPolygons);
    centerslength = length(centers);
    areas = zeros(centerslength, 1);
    triangleX = zeros(3, centerslength); 
    triangleY = zeros(3, centerslength);
    for i = 1:centerslength
    currentPoint = centers(i, :);
    distances = sqrt(sum((filteredPolygons - currentPoint).^2, 2));
    [~, closestIndices] = mink(distances, 3);
    triangleVertices = filteredPolygons(closestIndices, :);
    triangleX(:, i) = triangleVertices(:, 1); 
    triangleY(:, i) = triangleVertices(:, 2);
    triangle = polyshape(triangleVertices(:, 1), triangleVertices(:, 2));
    areas(i) = area(triangle);   
end

end

function [distances, penalty] = computeDistanceToVoronoi(cellX, cellY, particleX, particleY)
    % Ensuring points are inside the Voronoi cell by computing distance
    particleShape = polyshape(particleX, particleY);
    cellShape = polyshape (cellX, cellY);
    
    if isinterior(cellShape, particleX, particleY)
    penalty = 0; 
else
    % Checking if any vertices of the particle shape are outside the cell shape
    outsideVertices = ~isinterior(cellShape, particleX, particleY);
    
    if any(outsideVertices)
        penalty = 1e6; 
    end
    end
        distances = 0;
        
   
end

function err = objective(params, r, theta, cellX, cellY, initialX, initialY, targetArea, targetSolidFraction)
    dx = params(1); 
    dy = params(2);  
    scaleFactor = params(3);  
    

    scaled_rValues = scaleFactor * r(:,1);
    

    particleX = scaled_rValues .* cos(theta) + initialX + dx;
    particleY = scaled_rValues .* sin(theta) + initialY + dy;
    

    particleArea = pi * mean(scaled_rValues)^2;  
    
 
    [~, penalty] = computeDistanceToVoronoi(cellX, cellY, particleX, particleY);
    
  
    solidFractionError = (particleArea / targetArea) - targetSolidFraction;

    err = penalty + solidFractionError^2;
end


