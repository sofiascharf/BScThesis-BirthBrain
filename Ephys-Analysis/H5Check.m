
file = '2025-12-10T10-29-16McsRecording.h5';
h5disp(file);
datasetname = '/Data';

% Load
info = h5info(file, datasetname);
totalsamples = info.Datasets(1).Dataspace.Size;

% Read first min 
startRow = [1 1];
count = [60000 1];

lfp_segment = h5read(file, datasetname, startRow, count)