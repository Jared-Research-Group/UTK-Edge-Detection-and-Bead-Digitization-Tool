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

if sample_data.multibead 
    temp_size = size(sample_data.split);
    n_beads = temp_size(1);

    beadNum = 1;

    for i = 1:n_beads
        current_bead = sample_data.split{i};

        %Apply canny edge detection
        %edges =  edge(current_bead,'Canny',[0.35 0.75],25);
        edges =  edge(current_bead,'Canny',[0.4 0.65],35);

        fig1 = figure
        imshow(edges);
        beadName = strcat('Bead',num2str(beadNum));
        filepath = strcat(savepath,'\',beadName)
        saveas(fig1,strcat(filepath,' Edges.png'))
        beadNum = beadNum + 1;

        beadDigitization(edges,distperpix,filepath);
    end
else
    %Apply canny edge detection
    edges =  edge(analysis_sample,'Canny',[0.65 0.75],35);
    %edges =  edge(analysis_sample,'Canny',[0.4 0.65],35)

    figure
    imshow(edges);

    beadDigitization(edges,distperpix,filepath)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [] = beadDigitization(edges, distperpix,filepath)
%Get coordinates of boundary in user specified units
[imcoords.x, imcoords.y] = binary_coordinates(edges,distperpix);

%Plot digitized bead outline
fig2 = figure
fig2.Position(3:4) = [1900 400]
scatter(imcoords.x,imcoords.y)
title('Digitized Bead Outline')
xlabel('X distance (mm)')
ylabel('Y distance (mm)')

hold on

%Invert binary image 
I = double(~edges);

%Use AOFSkeletons Package to extract topological skeleton. 
[skeletonImage,fluxImage,distImage] = extract2DSkeletonFromBinaryImage(I);
%Call again to ensure skeleton is one to one line 
[skeletonImage,fluxImage,distImage] = extract2DSkeletonFromBinaryImage(skeletonImage);

%Get coordinates of binary topological skeleton 
[x_cent,y_cent] = binary_coordinates(skeletonImage,distperpix);
%Downsample the line to smooth out unwanted peaks and valleys
[x_cent,y_cent] = downsample(x_cent,y_cent,10);

%Use Pchip interpolation to create smooth line running through the
%downsampled topological skeleton 
xq = 0:0.1:max(imcoords.x);
pp = pchip(x_cent,y_cent,xq);
plot(xq,pp);

%Define endcap cutoff distance in mm 
start_cutoff = 4.5;
end_cutoff = max(xq)-4.5;
cutoffs = [start_cutoff,end_cutoff];

%Plot ideal centerline 
hold on
y_ctr = ideal_centerline(xq,pp);
plot(xq,y_ctr)

saveas(fig2,strcat(filepath,' Outline+Centerline.png'))

%Get width profile and plot
[width, start_idx, end_idx] = get_width(cutoffs,xq,pp,imcoords);
width = width.'
length = xq(start_idx:end_idx).';
fig3 = figure
plot(length,width)
xlabel('Length (mm)')
ylabel('Width (mm)')

saveas(fig3,strcat(filepath,' Width.png'));
T = table(length,width)
writetable(T, strcat(filepath, ' Width Profile.csv'))


end


function y_ctr = ideal_centerline(xq,pp)

y_ctr = pp(1)*ones(1,length(xq));


end

function std = ctr_line_deviation(pp,xq)


end

function [width, start_idx, end_idx] = get_width(cutoffs,xq,pp,imcoords)
[minstart, start_idx] = min(abs(xq - cutoffs(1)));
[minend, end_idx] = min(abs(xq - cutoffs(2)));

j = 1;

for i = start_idx:end_idx 
 
    diff = imcoords.x - xq(i);
    [x_min, first_idx] = min(abs(diff));
    idx = find(imcoords.x == imcoords.x(first_idx));
    
    for k = 1:length(idx)
     
        if imcoords.y(idx(k)) > pp(i) && length(idx) > 1
            y_high(j) = imcoords.y(idx(k));
        elseif imcoords.y(idx(k)) < pp(i) && length(idx) > 1
            y_low(j) = imcoords.y(idx(k));
        end
        
    end

    if length(idx) > 1
        j = j + 1;
    else
        xq(i) = [];
        end_idx = end_idx - 1;
        
    end
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
xmin = min(x);
x = x - xmin; 
end

function [x_new,y_new] = downsample(x,y,spacing)
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

