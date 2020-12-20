function [time,values] = readPitchData(fileIn)

% Author: Bart Doekemeijer, date: July 3, 2019 

%

% Description: this function reads the turbine blade pitch angles when

% defined as a three-dimensional vector (the IPC branch of SOWFA)

%



% Read the string data from the file

tmpTime{1} = [];

tmpData{1} = [];

file = fopen(fileIn,'r');

inputText = textscan(file,'%s',1,'delimiter','\n'); % Skip first line

    while (~feof(file))

        inputText = textscan(file,'%s',1,'delimiter','\n');

        inputText = inputText{1}{1};

        if length(inputText) > 2 % If line is not empty

            

            % Export turbine number and timestamp

            splitInput = strsplit(inputText,' ');

            turbineId = str2double(splitInput{1}) + 1;

            timeStamp = str2double(splitInput{2});

            

            % Export blade pitch data

            curlyBrcktIdx = strfind(inputText,'{');

            normalBrcktIdx = strfind(inputText,'(');

            if length(curlyBrcktIdx) > 0

                curlyBrcktEndIdx = strfind(inputText,'}');

                pitchString = inputText((curlyBrcktIdx+1):(curlyBrcktEndIdx-1));

                pitchValues = repmat(str2double(pitchString),1,3); % copy to each blade

            elseif length(normalBrcktIdx) > 0

                normalBrcktEndIdx = strfind(inputText,')');

                pitchString = inputText((normalBrcktIdx+1):(normalBrcktEndIdx-1));

                pitchValues = str2double(string(strsplit(pitchString,' ')));

            else
Q=str2num(inputText);
                pitchValues = [Q(4) Q(4) Q(4)]; % Write NaNs

            end

            

            % Write to a struct

            if turbineId > length(tmpTime)

                tmpTime{turbineId} = [];

                tmpData{turbineId} = [];

            end

            tmpTime{turbineId} = [tmpTime{turbineId}; timeStamp];

            tmpData{turbineId} = [tmpData{turbineId}; pitchValues];

        end

    end

    

% Convert to relevant turbine information

time = tmpTime{1};

values = tmpData;

nTurbs = length(tmpTime);

if any(unique([tmpTime{:}]) ~= time)

    error('There is some inconsistency in the timestamps.');

end

end