function padNumberedFiles(numberOfDigits,  nameForNewFiles, fileExtension, directory)
%Processing for various number of inputs
if nargin < 4
    directory = uigetdir(pwd, 'Select a folder');
end
if nargin < 3
    viableExtensions = {'.txt', '.pdf', '.stl', 'other'};
    [indx,tf] = listdlg('PromptString',{'Select File Extension you want to use'},...
    'SelectionMode','single','ListString',viableExtensions);
    if ~tf || strcmp('other', viableExtensions{indx})
        input = inputdlg('Please enter the file extension to use (eg: .txt)');
        fileExtension = input{1};
    else
        fileExtension = viableExtensions{indx};
    end
end
if nargin < 2
    yesNoAnswer = questdlg('Would you like to change the name of the file?');
    if(strcmp(yesNoAnswer, 'Yes'))
        input = inputdlg('Enter the file name: ');
        nameForNewFiles = input{1};
    else
        nameForNewFiles = '';
    end
end
if nargin < 1
    ranges = {'1-9', '10-99', '100 -999', '1000-9999', '10000+'};
    [indx,tf] = listdlg('PromptString',{'How many files to you have'},...
    'SelectionMode','single','ListString', ranges);
    if ~tf || indx >4
        input = inputdlg('How many digits would you like to pad to?\nWarning this may take a while with more than 10,000 files')
        numberOfDigits = str2double(input{1});
    else
        numberOfDigits = indx;
    end
end
% actual file renaming code
oldFolder = cd(directory);    
fileSearchString = strcat('*', fileExtension);
files = dir(fileSearchString);
padTextString = strcat('%0', num2str(numberOfDigits), 'd', fileExtension);
if isempty(nameForNewFiles)
    nameForNewFiles = files(1).name(1:end-numberOfDigits-length(fileExtension));
end
% Loop through each
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~, f] = fileparts(files(id).name);
    numberPart = f(end-numberOfDigits:end);
    %convert number string to actual number
    temp = numberPart -'0';
    index = find(temp >= 0 & temp < 10, 1);
    numberPart = str2double(numberPart(index:end));
    newFileName = strcat(nameForNewFiles, sprintf(padTextString, numberPart));
    % Convert to number
    % Pad the numbers in the name
    if ~strcmp(files(id).name, newFileName)
        movefile(files(id).name, newFileName);
    end
end
cd(oldFolder)
fprintf("The operation was performed successfully\n")
end
