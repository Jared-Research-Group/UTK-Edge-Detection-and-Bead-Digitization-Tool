%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bead Digitization Tool 
%Eduardo Miramontes, 2023
%emiramon@vols.utk.edu
%
%This program uses canny edge detection to extract the outline of a bead in
%binary format. A spatial calibration routine determines the distance per
%pixel based on a known distance entered by the user. This calibrated
%distance allows the bitmap of the bead outline to be plotted in a
%coordinate plane. The centerline of the outline is then computed. Once the
%outline and centerline are obtained, various kinds of analysis can be
%performed, such as generating a width profile, measurement of crookedness
%and waviness, etc. An auto-cropping tool that splits images of a plate 
% containing multiple beads into their respective beads is included, which
% enables batch analysis of individual beads. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear terminal window
clear
clc

% Read in a standard MATLAB gray scale demo image.
% Get the name of the file that the user wants to use.
defaultFileName = fullfile(cd, '*.*');
[baseFileName, folder] = uigetfile(defaultFileName, 'Select an image file');


% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
	% File doesn't exist -- didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end

%Read the image from file
analysis_sample = imread(fullFileName); 

%Turn analysis sample into grayscale image
analysis_sample =rgb2gray(analysis_sample);

%Apply spatial calibration
sample_data = spatial_calibration(analysis_sample);
distperpix = sample_data.distancePerPixel;

savepath = uigetdir(pwd,'Select folder to save results to')
beadNum = 26;
if sample_data.multibead 
    temp_size = size(sample_data.split);
    n_beads = temp_size(1);

    for i = 1:n_beads
        current_bead = sample_data.split{i};

        %Apply canny edge detection
        %edges =  edge(current_bead,'Canny',[0.35 0.75],25);
        edges =  edge(current_bead,'Canny',[0.4 0.5],35);

        fig1 = figure
        imshow(edges);
        beadName = strcat('Bead',num2str(beadNum));
        filepath = strcat(savepath,'\',beadName)
        saveas(fig1,strcat(filepath,' Edges.png'))
        beadDigitization(edges,distperpix,filepath,beadNum);
        beadNum = beadNum + 1;

    end
else
    %Apply canny edge detection
    %analysis_sample = imrotate(analysis_sample,-90)
    edges =  edge(analysis_sample,'Canny',[0.05 0.75],35);
    %edges =  edge(analysis_sample,'Canny',[0.4 0.65],35)
    imshow(analysis_sample);

    fig1 = figure
    imshow(edges);
    beadName = strcat('Bead',num2str(beadNum))
    filepath = strcat(savepath,'\',beadName)
    saveas(fig1,strcat(filepath,' Edges.png'))
    [imcoords] = beadDigitization(edges,distperpix,filepath,beadNum)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [imcoords] = beadDigitization(edges, distperpix,filepath,beadNum)
tablepath = 'C:\Users\emiramon\Documents\Data\Dataset 4\Start end points.csv'
startEndPts = readtable(tablepath);
startPts = startEndPts.Start;
endPts = startEndPts.End;

%Get coordinates of boundary in user specified units
[imcoords.x, imcoords.y] = binary_coordinates(edges,distperpix);

%Invert binary image 
I = double(~edges);

%Use AOFSkeletons Package to extract topological skeleton. 
[skeletonImage,fluxImage,distImage] = extract2DSkeletonFromBinaryImage(I);
%Call again to ensure skeleton is one to one line 
[skeletonImage,fluxImage,distImage] = extract2DSkeletonFromBinaryImage(skeletonImage);

%Get coordinates of binary topological skeleton 
skelPoints = [];
[skelPoints(1,:),skelPoints(2,:)] = binary_coordinates(skeletonImage,distperpix);


%Compute end point for center line on right side 
x_max = max(imcoords.x);
x_min = min(imcoords.x);
idx_max = find(imcoords.x == x_max(1));
mid_idx = round((idx_max(end) + idx_max(1))/2);
y_max = imcoords.y(mid_idx);

if isempty(find(skelPoints(1,:) == x_min))
    %Compute end point for center line on left side 
    idx_min = find(imcoords.x == x_min);
    mid_idx = round((idx_min(end) + idx_min(1))/2);
    y_min = imcoords.y(mid_idx);
end

%Downsample the line to smooth out unwanted peaks and valleys
skelPoints = downsample(skelPoints(1,:),skelPoints(2,:),3);
skelPoints(1,end) = x_max(1);
skelPoints(2,end) = y_max;
skelPoints(1,1) = x_min(1);
skelPoints(2,1) = y_min;

%Compute cutoffs for start and end points of torch
travelStart = startPts(beadNum);
travelEnd = endPts(beadNum);


%Use Pchip interpolation to create smooth line running through the
%downsampled topological skeleton 
centerLinePoints = [];
centerLinePoints(1,:) = travelStart:0.005:travelEnd;
centerLinePoints(2,:) = pchip(skelPoints(1,:),skelPoints(2,:),centerLinePoints(1,:));

%Compute ideal center line 
[y_ctr,angle] = ideal_centerline(centerLinePoints);
idealCtrLinePts = [centerLinePoints(1,:);y_ctr];
disp(size(idealCtrLinePts))

%Calculate transformation matrix and apply transformation to rotate
%misaligned beads so that ideal center line is parallel to x axis
impoints = [imcoords.x.';imcoords.y.'];
T = [cosd(angle),sind(angle);-sind(angle), cosd(angle)];
impoints = T*impoints;
imcoords.x = impoints(1,:);
imcoords.y = impoints(2,:);
clear impoints

centerLinePoints = T*centerLinePoints;
idealCtrLinePts = T*idealCtrLinePts;
skelPoints = T*skelPoints;
startPts = T(1,1)*startPts;
endPts = T(1,1)*endPts;

%Shift points so that bead and centerline start at 0 
% xmin = min(imcoords.x);
% imcoords.x = imcoords.x - xmin;
% startPts = startPts - xmin;
% endPts = endPts - xmin;
% skelPoints(1,:) = skelPoints(1,:) - xmin;
% centerLinePoints(1,:) = centerLinePoints(1,:) - xmin;
% idealCtrLinePts(1,:) = idealCtrLinePts(1,:) - xmin;


%Plot digitized bead outline
fig2 = figure
fig2.Position(3:4) = [4000 300];
scatter(imcoords.x,imcoords.y)
title('Digitized Bead Outline')
xlabel('X distance (mm)')
ylabel('Y distance (mm)')
hold on 
scatter(skelPoints(1,:),skelPoints(2,:));
hold on 
plot(centerLinePoints(1,:),centerLinePoints(2,:));
saveas(fig2,strcat(filepath,' Outline.png'))

fig3 = figure
set(gcf,'Position',[100 100 4000 300])
%Plot center line 
%scatter(skelPoints(1,:),skelPoints(2,:));
%hold on
scatter(centerLinePoints(1,:),centerLinePoints(2,:));
%Plot ideal centerline  
hold on
plot(idealCtrLinePts(1,:),idealCtrLinePts(2,:))
xlabel('X (mm)')
ylabel('Y (mm)')
saveas(fig3,strcat(filepath,' Centerline.png'))

%Get width profile 
cutoffs(1) = startPts(beadNum);
cutoffs(2) = endPts(beadNum);
Width = get_width(cutoffs,centerLinePoints,imcoords);
Width = Width.';
Length = centerLinePoints(1,:).';

%Shift length so it starts at 0 
minLength = min(Length)
Length = Length - minLength;

Length = Length(10:end);
Width = Width(10:end);

%Plot width
fig4 = figure
plot(Length,Width)
xlabel('Length (mm)')
ylabel('Width (mm)')
fig4.Position(3:4) = [4000 300];

%Compute deviation of actual center line from ideal center line 
CenterLineDeviation = ctr_line_deviation(centerLinePoints,idealCtrLinePts);

CenterLineDeviation = CenterLineDeviation(10:end).';

saveas(fig4,strcat(filepath,' Width.png'));
T = table(Length,Width,CenterLineDeviation);
writetable(T, strcat(filepath, ' Width Profile.csv'))


end


function [y_ctr, angle] = ideal_centerline(centerLinePoints)

xq = centerLinePoints(1,:);
pp = centerLinePoints(2,:);

m = (pp(end) - pp(1))/(xq(end) - xq(1));
y_ctr = m * (xq - xq(1)) + pp(1);

angle = atand(m);

end

function difference = ctr_line_deviation(centerLinePoints,idealCtrLinePts)

difference = centerLinePoints(2,:) - idealCtrLinePts(2,:);

end

function width = get_width(cutoffs,centerLinePoints,imcoords)

xq = centerLinePoints(1,:);
pp = centerLinePoints(2,:);

% [minstart, start_idx] = min(abs(xq - cutoffs(1)));
% [minend, end_idx] = min(abs(xq - cutoffs(2)));

j = 1;

for i = 1:length(xq)

diff = imcoords.x - xq(i);
[sortedDiff,sortedIdx] = sort(abs(diff));
xmin = sortedDiff(1);
first_idx = sortedIdx(1);
idx = find(imcoords.x == imcoords.x(first_idx));

%Boolean that stores info on whether current index is above or below
%centerline. Above is true, below is false 

if length(idx) > 2

    yvals = imcoords.y(idx);
    yhighsidx = find(yvals > pp(i));
    ylowsidx = find(yvals < pp(i));

    yhighs = yvals(yhighsidx)
    newyhigh = min(yhighs)
    newhighidx = find(yvals == newyhigh)

    ylows = yvals(ylowsidx)
    newylow = max(ylows)
    newlowidx = find(yvals == newylow)
    
    idx = [idx(newhighidx) idx(newlowidx)]

end

for n = 1:length(idx)

    direction = false;
    if imcoords.y(idx(n)) > pp(i)
        y_high(j) = imcoords.y(idx(n));
        direction = true;
    elseif imcoords.y(idx(n)) < pp(i)
        y_low(j) = imcoords.y(idx(n));
    end
    
    found = false;

    if length(idx) > 1
        found = true;
    end

    newDirection = false;
    k = 2;

    while (found == false)

        second_idx = sortedIdx(k);
        testValY = imcoords.y(second_idx);

        if (testValY > pp(i))
            newDirection = true;
        else
            newDirection = false;
        end

        if (direction ~= newDirection)
            found = true;
        end

        k = k + 1;

    end

    if newDirection && length(idx) == 1 
        y_high(j) = testValY;
        j = j + 1;
    elseif newDirection == false && length(idx) == 1
        y_low(j) = testValY;
        j = j + 1;
    elseif n == 2
        j = j + 1;
    end

end


%     for k = 1:length(idx)
%      
%         if imcoords.y(idx(k)) > pp(i) && length(idx) > 1
%             y_high(j) = imcoords.y(idx(k));
%         elseif imcoords.y(idx(k)) < pp(i) && length(idx) > 1
%             y_low(j) = imcoords.y(idx(k));
%         end        
%     end
% 
%     if length(idx) > 1
%         j = j + 1;
%     else
%         xq(i) = [];
%         end_idx = end_idx - 1;
%         
%     end
end
width = y_high - y_low;
end


function [x, y] = binary_coordinates(binaryImage,distperpix)
%Function applies distance per pixel to indexes on the binary image that
%contain the boundary. Each matrix element represents a pixel. Row indexes
%correspond to x coordinates and column indexes correspond to y coordinates

%Find indexes on binary image that contain a 1
[row,col] = find(binaryImage);

%Multiply indexes by distance per pixel 
y = row*distperpix;
x = col*distperpix;

%Shift all x values by the minimum x to make plot start at x = 0 
end

function newSkel = downsample(x,y,spacing)
%Downsample
x_new = [];
y_new = [];
j = 1;
for i = 1:length(x)-1
    if mod(i,spacing) == 0 && i >= spacing
        x_new(j) = x(i);
        y_new(j) = y(i);
        j = j + 1;
    end
end
x_new(end) = x(end);
y_new(end) = y(end);

newSkel = [];
newSkel(1,:) = x_new;
newSkel(2,:) = y_new;

end

function [x_new,y_new] = thin_centerline(x,y)

j = 1; 

for i = 2:length(x)
    if x(i) == x(i-1) 
        x_new(j) = x(i);
        y_new(j) = (y(i) + y(i-1))/2;
        i = i + 1
        j = j + 1;
    elseif  i+1 <= length(x)
        if x(i) ~= x(i+1)
            x_new(j) = x(i);
            y_new(j) = y(i);
            j = j + 1;  
            i = i + 1;
        end

    end
end


end

function [] = plot_bead_centerline(vx_sorted,vy_sorted,imcoords)
figure 
scatter(imcoords.x,imcoords.y)
hold on 
scatter(vx_sorted,vy_sorted)
xlabel('X distance (mm)')
ylabel('Y distance (mm)')

end

