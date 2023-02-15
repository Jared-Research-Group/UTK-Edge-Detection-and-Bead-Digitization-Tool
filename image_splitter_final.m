 % This script splits images with multiple beads into seperate images based on the number of beads, width of each bead, and the length of the image. It assumes that the beads are equally spaced! It operates independently from 'spatial_calibration_main' and 'spatial_calibration', so that it can be used to split multiple images and simply load them into memory.
% !!! This script does not automatically load the newly split images into any other script !!!

% Adry Lain 12/20/2022
% alain@vols.utk.edu

% Clear the terminal window.
clc;
% Open the workspace so the user can visually confirm the image was split and loaded.
workspace;
fontSize = 20;
success = false;
yoffset = 0;

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

% Read the image from file.
multibead = imread(fullFileName);
% Get the dims of the image.
[rows, columns, numberOfColorBands] = size(multibead);

%Main while loop.
while success == false 

% Initial image display.

figureHandle = figure;
subplot(1,2,1);
imshow(multibead, []);
axis on;
title('Original Multi-Bead Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);


% Ask user to start.
message = sprintf('Begin splitting multi-bead image?');
reply = questdlg(message, 'Split', 'OK', 'Cancel', 'OK');
if strcmpi(reply, 'Cancel')
	% User said Cancel, so exit.
	return;
end

% Ask about orientation of the image.

message = sprintf('Is the image oriented with the bead length going side to side?');
reply = questdlg(message, 'Orientation', 'Yes', 'No', 'Yes');

if strcmpi(reply, 'No')
	% For ease of use, the program will now auto rotate the image.
	multibead = imrotate(multibead, 90);
	% Show new image.
	imshow(multibead, []);
	axis on;
	title('Reoriented Multi-Bead Image', 'FontSize', fontSize);
end


% Prompt for number of beads.
userPrompt = {'Enter the number of individual beads in the image:'};
dialogTitle = 'Specify number of beads';
numberOfLines = 1;
def = {'beads'};
numberOfBeads = inputdlg(userPrompt, dialogTitle, numberOfLines, def);
if isempty(numberOfBeads)
	return;
end
numberOfBeads = str2double(numberOfBeads);


% Prompt to draw a line to record bead width.

message = sprintf('Draw a line from one top edge of any bead to its bottom edge in order to record its width. \nFirst, left-click to anchor the first point.\nRight-click or double-left-click to anchor the second endpoint of the line.'); 
reply = questdlg(message, 'Bead Width', "OK", 'Cancel', 'OK');

if strcmpi(reply, 'Cancel')
	% User would like to cancel.
	return;
end


[beadVerticalDistance,yoffset] = MeasureVerticalPixelDistance();

% Prompt for start and end position line for the entire image width.

message = sprintf('Now draw a line from the top edge of the first bead to the bottom edge of the last bead.'); 
reply = questdlg(message, 'Boundaries', "OK", 'Cancel', 'OK');

if strcmpi(reply, 'Cancel')
	% User would like to cancel.
	return;
end

[multibeadVerticalDistance,yoffset] = MeasureVerticalPixelDistance();

% Use the information to cut up the multibead image and store the pieces in seperate tables.

% Calculate the total amount of empty space between beads.
totalMarginSpace = multibeadVerticalDistance - (beadVerticalDistance*(numberOfBeads));
% Calculate what the total margin space for each individual image should be.
marginSpace = totalMarginSpace/(numberOfBeads-1);
% Calculate how tall each individual image should be.
individualImageHeight = beadVerticalDistance + marginSpace;
% Cut em' up! 
splitImages = cell(numberOfBeads,1);
for i = 1:numberOfBeads

	% Add yoffset to the product of the index by the image height.
	bottomLocation = floor((yoffset+i*individualImageHeight));
	% Same thing, but subtract one image height to get the top.
	topLocation = floor((yoffset+(i*individualImageHeight-individualImageHeight)));

	splitImages{i} = multibead(topLocation:bottomLocation, 1:end, 1:end);
end

% Display an image sample on the screen, and ask the user if it looks correct; otherwise, restart the process.

subplot(1,3,1);
sampleImage = cell2mat(splitImages(1,1));
imshow(sampleImage, []);
axis on;
title('Sample Result', 'FontSize', fontSize);

message = sprintf('Does this sample look correct?'); 
reply = questdlg(message, 'Sample Result', "Yes", 'No', 'Yes');

if strcmpi(reply, 'Yes')
	% The user is satisfied with the result. End loop.
	success = true;
end

end % End of main while loop.


function offset = setYoffset()

end


% Rewritten and seperate from the code in "spatial_calibration" due to a need to make a class in order to call it and un-necessary code (for the purposes of this script).
function [distanceInPixels,yoffset] = MeasureVerticalPixelDistance()

	subplot(1,2, 1); % Switch to image axes.
	[interpolatedXCoords,interpolatedYCoords, rgbValues, userClicked_xi,userClicked_yi] = improfile(1000);
	% Quit if they didn't click at least two points.
	if length(userClicked_xi) < 2
		return;
	end
		
	% rgbValues is 1000x1x3.  Call Squeeze to get rid of the singleton dimension and make it 1000x3.
	rgbValues = squeeze(rgbValues);
	% Get the VERTICAL distance from the first clicked point to the second clicked point, in pixels.
	distanceInPixels = abs(userClicked_yi(2)-userClicked_yi(1));

	% Grab the highest y position to use as the offset when splitting. Since the total image height operation is called last, the data recorded for the bead vertical distance operation will be overwritten and not used.
	if userClicked_yi(2) < userClicked_yi(1)
		yoffset = userClicked_yi(2);
	else
		yoffset = userClicked_yi(1);
	end

	% Plot the line visually.
	hold on;
	lastDrawnHandle = plot(userClicked_xi, userClicked_yi, 'y-', 'LineWidth', 2);
	

end
