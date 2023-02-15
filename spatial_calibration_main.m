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
bead = imread(fullFileName); 

%Turn bead into grayscale image 
bead=rgb2gray(bead);

%Apply spatial calibration
calibration = spatial_calibration(bead);
distperpix = calibration.distancePerPixel;

