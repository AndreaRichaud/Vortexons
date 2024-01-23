function [] = Produce_animation(filename)
% Specify the folder where the figures are stored
folderPath = 'Temp_real';

% Specify the output video file name
videoFileName = sprintf('%s.avi',filename);

% Create a VideoWriter object
videoWriter = VideoWriter(videoFileName);
videoWriter.FrameRate = 15; % Set the frame rate (adjust as needed)

% Open the video writer
open(videoWriter);

% Get the list of figure files in the folder
figureFiles = dir(fullfile(folderPath, '*.jpg')); % Change '*.png' to the actual file format of your figures

% Loop through each figure file and add it to the video
for i = 1:length(figureFiles)
    % Read the figure
    currentFigure = imread(fullfile(folderPath, figureFiles(i).name));
    
    % Add the frame to the video
    writeVideo(videoWriter, currentFigure);
end

% Close the video writer
close(videoWriter);

end

