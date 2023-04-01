% This demo allows you to spatially calibrate your image and then make distance or area measurements.

function output = spatial_calibration(originalImage)

% Check that user has the Image Processing Toolbox installed.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

hasIPT = license('test', 'image_toolbox');
if ~hasIPT
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
end

% Get the dimensions of the image.
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(originalImage);
% Display the original gray scale image.
figureHandle = figure;
subplot(1,2, 1);
imshow(originalImage, []);
axis on;
title('Original Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf,'name','Demo by ImageAnalyst','numbertitle','off')

button = 99;

while button ~= 6
	if button > 1
		% Let them choose the task, once they have calibrated.
		button = menu('Select a task', 'Manual Calibration','External Calibration', 'Measure Distance','Split Beads','Speed Override','Exit');
	end
	switch button
		case 1
			[success,output] = Calibrate();
			% Keep trying if they didn't click properly.
			while ~success
				[success,output] = Calibrate();
			end
			% If they get to here, they clicked properly
			% Change to something else so it will ask them
			% for the task on the next time through the loop.
			button = 99;
		case 2
            output = ExternalCalibration();
			
        case 3
            DrawLine();            

        case 4

            output = SplitBeads(originalImage, output);

        case 5
            output.units = 'mm';
	        output.distancePerPixel = 0.047256;
            output.multibead = false
        otherwise
           
			close(figureHandle);
			break;
	end
end

end

%=====================================================================
function [success,calibration] = Calibrate()
try
	success = false;
    
	[cx, cy, rgbValues, xi,yi] = improfile(1000);
	% rgbValues is 1000x1x3.  Call Squeeze to get rid of the singleton dimension and make it 1000x3.
	rgbValues = squeeze(rgbValues);
	distanceInPixels = sqrt( (xi(2)-xi(1)).^2 + (yi(2)-yi(1)).^2);
	if length(xi) < 2
		return;
	end
	% Plot the line.
	hold on;
	lastDrawnHandle = plot(xi, yi, 'y-', 'LineWidth', 2);
	
	% Ask the user for the real-world distance.
	userPrompt = {'Enter real world units (e.g. microns):','Enter distance in those units:'};
	dialogTitle = 'Specify calibration information';
	numberOfLines = 1;
	def = {'mm', '500'};
	answer = inputdlg(userPrompt, dialogTitle, numberOfLines, def);
	if isempty(answer)
		return;
	end
	calibration.units = answer{1};
	calibration.distanceInPixels = distanceInPixels;
	calibration.distanceInUnits = str2double(answer{2});
	calibration.distancePerPixel = calibration.distanceInUnits / distanceInPixels;
	success = true;
    calibration.multibead = false
	
	message = sprintf('The distance you drew is %.2f pixels = %f %s.\nThe number of %s per pixel is %f.\nThe number of pixels per %s is %f',...
		distanceInPixels, calibration.distanceInUnits, calibration.units, ...
		calibration.units, calibration.distancePerPixel, ...
		calibration.units, 1/calibration.distancePerPixel);
	uiwait(msgbox(message));
catch ME
	errorMessage = sprintf('Error in function Calibrate().\nDid you first left click and then right click?\n\nError Message:\n%s', ME.message);
	fprintf(1, '%s\n', errorMessage);
	WarnUser(errorMessage);
end


return;	% from Calibrate()
end

function calibration = ExternalCalibration()

% Ask the user for the pixel calibration from external source 
userPrompt = {'Enter pixel calibration from external source (units/pixel): ',
    'Enter distance units: '};
dialogTitle = 'Specify calibration information';
numberOfLines = 1;
def = {'0.05', 'mm'};
answer = inputdlg(userPrompt, dialogTitle, numberOfLines, def);

calibration.units = answer{2};
calibration.distancePerPixel = str2double(answer{1})
calibration.multibead = false

return;

end


%=====================================================================
% --- Executes on button press in DrawLine.
function DrawLine()
try
	fontSize = 14;
	
	instructions = sprintf('Draw a line.\nFirst, left-click to anchor first endpoint of line.\nRight-click or double-left-click to anchor second endpoint of line.\n\nAfter that I will ask for the real-world distance of the line.');
	title(instructions);
	msgboxw(instructions);
	subplot(1,2, 1); % Switch to image axes.
	[interpolatedXCoords,interpolatedYCoords, rgbValues, userClicked_xi,userClicked_yi] = improfile(1000);
	% Quit if they didn't click at least two points.
	if length(userClicked_xi) < 2
		return;
	end
	
	% Get the profile again but spaced at the number of pixels instead of 1000 samples.
	theImage = getimage(gca); % Extract image out of axes control.
	% Find the distance in pixels from point 1 to point 2.
	lineLength = round(sqrt((userClicked_xi(1)-userClicked_xi(2))^2 + (userClicked_yi(1)-userClicked_yi(2))^2))
	[interpolatedXCoords,interpolatedYCoords, rgbValues] = improfile(theImage, userClicked_xi, userClicked_yi, lineLength);
	
	% rgbValues is 1000x1x3.  Call Squeeze to get rid of the singleton dimension and make it 1000x3.
	rgbValues = squeeze(rgbValues);
	% Get the distance from the first clicked point to the second clicked point, in pixels.
	distanceInPixels = sqrt( (userClicked_xi(2)-userClicked_xi(1)).^2 + (userClicked_yi(2)-userClicked_yi(1)).^2);
	% Now convert that to real world units.
	distanceInRealUnits = distanceInPixels * calibration.distancePerPixel;
	
	% Plot the line.
	hold on;
	lastDrawnHandle = plot(userClicked_xi, userClicked_yi, 'y-', 'LineWidth', 2);
	
	% Plot profiles along the line of the red, green, and blue components.
	subplot(1,2,2);
	[rows, columns] = size(rgbValues);
	if columns == 3
		% It's an RGB image.
		plot(rgbValues(:, 1), 'r-', 'LineWidth', 2);
		hold on;
		plot(rgbValues(:, 2), 'g-', 'LineWidth', 2);
		plot(rgbValues(:, 3), 'b-', 'LineWidth', 2);
		title('Red, Green, and Blue Profiles along the line you just drew.', 'FontSize', 14);
	else
		% It's a gray scale image.
		plot(rgbValues, 'k-', 'LineWidth', 2);
	end
	xlabel('X', 'FontSize', fontSize);
	ylabel('Gray Level', 'FontSize', fontSize);
	title('Intensity Profile', 'FontSize', fontSize);
	grid on;
	
	% Inform user via a dialog box.
	txtInfo = sprintf('Distance = %.1f %s, which = %.1f pixels.', ...
		distanceInRealUnits, calibration.units, distanceInPixels);
	msgboxw(txtInfo);
	% Print the values out to the command window.
	fprintf(1, '%\n', txtInfo);
	
catch ME
	errorMessage = sprintf('Error in function DrawLine().\n\nError Message:\n%s', ME.message);
	fprintf(1, '%s\n', errorMessage);
	WarnUser(errorMessage);
end
end  % of DrawLine()

function output = SplitBeads(multibead,output)

workspace;
fontSize = 20;
success = false;
yoffset = 0;
debugCircleSize = 5;

%Main while loop.
while success == false 

topTruncate = false;
bottomTruncate = false;
debug = false;
quick = false;
newUser = false;

% Ask for mode.

% message = sprintf('Select New User if you do not know how to use this tool.');
% reply = questdlg(message, 'Mode', 'Quick', 'New User', 'Debug', 'Quick');
% 
% switch reply
% 	case 'Quick'
% 		quick = true;
% 	case 'New User'
% 		newUser = true;
% 	case 'Debug'
% 		debug = true;
% end

%Hard coded for speed during debugging
quick = true 

% Initial image display.

imshow(multibead, []);
axis on;
title('Original Multi-Bead Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Ask about orientation of the image.
% message = sprintf('Is the image oriented with the bead length going left to right (horizontal)?');
% reply = questdlg(message, 'Orientation', 'Yes', 'No', 'Yes');
% 
% if strcmpi(reply, 'No')
% 	% For ease of use, the program will now auto rotate the image.
% 	multibead = imrotate(multibead, -90);
%     hold off
% 	% Show new image.
% 	imshow(multibead, []);
% 	axis on;
% 	title('Reoriented Multi-Bead Image', 'FontSize', fontSize);
% end

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

if quick == false
message = sprintf('Draw a line from one top edge of any bead to its bottom edge in order to record its width. \nFirst, left-click to anchor the first point.\nRight-click or double-left-click to anchor the second endpoint of the line.'); 
reply = questdlg(message, 'Bead Width', "OK", 'Cancel', 'OK');

if strcmpi(reply, 'Cancel')
	% User would like to cancel.
	return;
end
end

[beadVerticalDistance,yoffset] = MeasureVerticalPixelDistance();

% Prompt for start and end position line for the entire image width.

if quick == false
message = sprintf('Now draw a line from the top edge of the first bead to the bottom edge of the last bead.'); 
reply = questdlg(message, 'Boundaries', "OK", 'Cancel', 'OK');

if strcmpi(reply, 'Cancel')
	% User would like to cancel.
	return;
end
end

[multibeadVerticalDistance,yoffset] = MeasureVerticalPixelDistance();

% Use the information to cut up the multibead image and store the pieces in seperate tables.

% Calculate the total amount of empty space between beads.
totalMarginSpace = multibeadVerticalDistance - (beadVerticalDistance*(numberOfBeads));
% Calculate what the total margin space for each individual image should be.
marginSpace = totalMarginSpace/numberOfBeads;
% Calculate how tall each individual image should be.
individualImageHeight = beadVerticalDistance + marginSpace;


% Find out if the top or the bottom boundary does not fit with the calculated individual image height.
if (yoffset < ((individualImageHeight-beadVerticalDistance)/2)) 
	topTruncate = true;
	fprintf("Needs to be truncated on the top!")
elseif rows < ((yoffset + numberOfBeads*individualImageHeight))
	bottomTruncate = true;	
	yoffset + numberOfBeads*individualImageHeight
	fprintf("Needs to be truncated on the bottom!")
end


% Cut em' up! 
splitImages = cell(numberOfBeads,1);
for i = 1:numberOfBeads

	if topTruncate == true && i == 1
	
	    topLocation = 1;
	    bottomLocation = floor((yoffset+i*individualImageHeight));
        splitImages{i} = multibead(topLocation:bottomLocation, 1:end, 1:end);
		
		if debug == true 
			plot(columns/2, topLocation, 'ro', 'MarkerSize', debugCircleSize);
			plot(columns/2, bottomLocation, 'ro', 'MarkerSize', debugCircleSize);
		end

	elseif bottomTruncate == true && i == numberOfBeads

	    topLocation = floor((yoffset+((i-1)*individualImageHeight)));
	    bottomLocation = rows;
        splitImages{i} = multibead(topLocation:bottomLocation, 1:end, 1:end);

		if debug == true 
			plot(columns/2, topLocation, 'ro', 'MarkerSize', debugCircleSize);
			plot(columns/2, bottomLocation, 'ro', 'MarkerSize', debugCircleSize);
		end
		
    else
	    % Add yoffset to the product of the index by the image height.
	    bottomLocation = floor((yoffset+i*individualImageHeight));
	    % Same thing, but subtract one image height to get the top.
	    topLocation = floor((yoffset+((i-1)*individualImageHeight)));
    
	    splitImages{i} = multibead(topLocation:bottomLocation, 1:end, 1:end);

		if debug == true 
			plot(columns/2, topLocation, 'ro', 'MarkerSize', debugCircleSize);
			plot(columns/2, bottomLocation, 'ro', 'MarkerSize', debugCircleSize);
		end
	
	end
end

output.split = splitImages;
output.multibead = true
%Pause after splitting images for visual markers.
if debug == true
	message = sprintf('DEBUG WAIT.');
	reply = questdlg(message, 'DEBUG WAIT', 'Proceed', 'Cancel', 'Proceed');
	if reply == 'Cancel'
		return;
	end
end


% Display an image sample on the screen, and ask the user if it looks correct; otherwise, restart the process.

subplot(1,3,1);
sampleImage = cell2mat(splitImages(1,1));
imshow(sampleImage, []);
axis on;
title('Sample Result', 'FontSize', fontSize);

% message = sprintf('Does this sample look correct? (This is the first image, if the image was truncated, the top may be cut short. This is normal.)'); 
% reply = questdlg(message, 'Sample Result', "Yes", 'No', 'Yes');

message = sprintf('Are there any samples you want to omit from the dataset?')
reply = questdlg(message,'Omit Samples?','Yes','No','Yes')

if strcmpi(reply,'Yes')

    userPrompt = {'Enter the index of the bead you want to omit (start counting with 1 from the top)'};
    dialogTitle = 'Specify beads to omit';
    numberOfLines = 1;
    def = {'Index'};
    omit_idx = inputdlg(userPrompt,dialogTitle,numberOfLines,def)
    omit_idx = str2num(omit_idx{1})

    for i = 1:length(omit_idx)
        idx = omit_idx(i)
        output.split{idx} = [];
    end

     output.split = output.split(~cellfun('isempty',output.split))
    
end

% message = sprintf('Are you satisfied with the results?')
% reply = questdlg(message, 'Exit Splitter', "Yes", 'No', 'Yes');
% 
% if strcmpi(reply, 'Yes')
% 	% The user is satisfied with the result. End loop.
% 	success = true;
% end

success = true;

debug = false;
quick = false;
newUser = false;

end % End of main while loop.

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
%=====================================================================
function msgboxw(message)
uiwait(msgbox(message));
end
%=====================================================================
function WarnUser(message)
uiwait(msgbox(message));
end