function getBDataEstudiantes(file_name)

data_folder = 'data'; % Assuming your script is in the parent directory of 'data'
data_file_folders = dir(fullfile(data_folder, 'datos*'));    
for i = 1:length(data_file_folders)
  file_path = fullfile(data_folder, data_file_folders(i).name, file_name);
    
% Get the file extension
[~, ~, ext] = fileparts(file_name);
%Create a conditional to read excel files
if strcmp(ext, '.xlsx') 
    %Read Excel
    T = readtable(file_path, 'Sheet', 1);
elseif strcmp(ext, '.csv')
    T = readtable(file_path);
%else
 %   error('Unsupported file format. Only .xlsx and .csv files are supported.');
end
% Extract variables of interest for analysis: subject#, session#,
% ss_am, ss_del, ll_am, ll_del, respT_corr, respT_rt.
% for respT.corr, 1=chose delayed and 0=chose immediate

data = [T.ss_am T.ss_del T.ll_am T.ll_del T.respT_corr T.respT_rt]; %output psychopy csv

% Remove practice trials
%data = data(10:end,:);g

% Remove rows with NaNs:
data(isnan(data(:,1)),:) = [];

% % Total number of trials imported - this number can change accordingly when
% % non-answered trials are removed below
trialNum = length(data);
sub = T.participant(1);
ses = T.session(1);

if iscell(sub)==1
    sub = str2double(sub);
end
if iscell(ses)==1
    ses = str2double(ses);
end

% if isnan(sub) = 1
%     sub = raw(10,35);
% end
% 
% if isnan(ses) = 1
%     ses = raw(10,34);
% end

thisdir = pwd;
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'bData'))
    disp('Making data directory');
    mkdir('bData');
end

% save output file as structure containing all fields of 'data' + subject
% and session number
outfile = fullfile('bData',sprintf('BE_S%d_%d.mat',sub,ses));
save(outfile,'data','sub','ses','trialNum');

end