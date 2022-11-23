%test_SOFAdbURL - Test script, using all SOFA files from SOFAdbURL database.
% 
%   test_SOFAdbURL scans all [SOFAdbURL '/database/'] subfolders for *.sofa files.
%   The SOFA files are downloaded to the local path [SOFAdbPath '\urlDatabase'], loaded to a SOFA object, and saved to a temporary file.
%   All warnings and error messages are stored in a log file in the local folder (log.csv).
%   Script only working in Matlab, not in Octave

% #Author: Michael Mihocic: header documentation updated (13.10.2022)
% #Author: Michael Mihocic: some updated, mostly regarding convention upgrades and stability (18.11.2022)
% #Author: Michael Mihocic: bugs fixed for sofa files creating more than one warning (22.11.2022)
% #Author: Michael Mihocic: stability improved (webread); file naming changed (include date, time) (23.11.2022)
% 
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% Clean up and set parameters
clear;
downloadAttempts = 8; % how often try to download until file is skipped
upgradeConventions = false; % true; % upgrade SingleRoomDRIR2SingleRoomSRIR and MultiSpeakerBRIR2SingleRoomMIMOSRIR?

%% Preparations
clc; close all; % clean-up first
tic; % timer
SOFAstart('restart');
% warning('query','all');
warning ('on','all');
warning('off','SOFA:save:API');

%% Prologue
disp('      ');
disp('############################################');
disp('########   TEST SOFA URL DATABASES  ########');
disp('############################################');
disp('      ');

%% log files
starttime=datestr(now(), 'yymmdd HHMMSS');
logfile=[SOFAdbPath '\urlDatabase\log ' starttime '.csv'];
logfileErr=[SOFAdbPath '\urlDatabase\logErr ' starttime '.csv'];
logIndex=1; logErrIndex=1; % start indices in log file
% need to create sub directory first?
if ~exist([SOFAdbPath '\urlDatabase\'], 'dir')
   mkdir([SOFAdbPath '\urlDatabase\']);
end
% log
fid=fopen(logfile,'w');
fprintf(fid, '%s\n\n%s',[ '*** Checking SOFA files for errors and warnings while downloading, loading & saving; start time: ' datestr(now(), 'dd.mm.yyyy - HH:MM:SS')],['INDEX' char(9) 'DB Index' char(9) 'TYPE' char(9) 'Operation' char(9) 'Message' char(9) 'File/Link']);
fclose(fid);
% log errors
fid=fopen(logfileErr,'w');
fprintf(fid, '%s\n\n%s',[ '*** Checking SOFA files for errors while downloading, loading & saving; start time: ' datestr(now(), 'dd.mm.yyyy - HH:MM:SS')],['INDEX' char(9) 'DB Index' char(9) 'TYPE' char(9) 'Operation' char(9) 'Message' char(9) 'File/Link']);
fclose(fid);
lastwarn(''); % clear last warning
errorCatch=false; % use variable to detect catched errors
disp(['### See log file for errors and warnings: ' logfile])
disp(['### See log file for errors: ' logfileErr])

%% Get subdirectories from URL db
dirDatabase=[SOFAdbURL '/database/'];
options = weboptions("Timeout", 60, "UserAgent", "Mozilla"); % Matlab
dirContent=webread(dirDatabase,options);

% Get subdirectories in database
subDirs = extractBetween(dirContent, "alt=""[DIR]""></td><td><a href=""", "/""");

% temporary target save file, overwritten over in every loop
temptargetfile = [SOFAdbPath '\urlDatabase\temp.sofa'];

%% Loop to scan, download and handle data
for ii=1:size(subDirs,1)
%     errorCatch=false; % ok, give it a chance

    % let's scan the current subdirectory
    currentSubLink=[dirDatabase char(subDirs(ii))];

    %% webread in directory (multiple attempts)
    attempts = 0; % count download attempts
    for retry=1:downloadAttempts
        try % to download the file
            errorCatch=false;
%             websave(targetfile,downloadLink,options);
            currentContent=webread(currentSubLink,options);
        catch me
            errorCatch=true;

            if retry < downloadAttempts % give it another try
                disp(['Webread attempt failed, retry webread (' num2str(retry) '/' num2str(downloadAttempts-1) ') from: ' currentSubLink]);
                appendToFile(logfile,[num2str(logIndex) char(9) num2str(ii) char(9) 'RETRY' char(9) 'Webread attempt failed for the ' num2str(retry) '. time, wait and retry...' char(9) currentSubLink]);   logIndex = logIndex + 1;
                pause(retry * 5); % give it a break of a couple of seconds, server might block too many simulataneous downloads!?
            else % ok, let's skip this file
%                 warning(['Download' char(9) me.message]);
                errorCatch=true;
                warning(['webread' char(9) me.message]);
                pause(5); % give it a break, server might block multiple simultaneous downloads!?
            end

        end
        if errorCatch==false
            break % leave loop
        end
    end

%     try
%         currentContent=webread(currentSubLink,options);
%     catch me
%         errorCatch=true;
%         warning(['webread' char(9) me.message]);
%         pause(5); % give it a break, server might block multiple simultaneous downloads!?
%     end

    if errorCatch==false
        % let's look for files in this subdirectory
        currentFiles = extractBetween(currentContent, "</td><td><a href=""", """>");
        % where do we wanna store the data?
        targetDir=[SOFAdbPath '\urlDatabase\' char(strrep(subDirs(ii),"%20", " "))];
    
        % a short status message:
        disp('      ');
        disp(['### Downloading and checking URL database ' num2str(ii) ' from ' num2str(size(subDirs,1)) ': ' char(strrep(subDirs(ii),"%20", " "))]);
        disp(['    from: ' currentSubLink]);
        disp(['    to:   ' targetDir]);
    
        % need to create sub directory first?
        if ~exist(targetDir, 'dir')
           mkdir(targetDir);
        end
        
        % loop for every file
        for jj=1:size(currentFiles,1)
            downloadLink=[currentSubLink '/' char(currentFiles(jj))];
            if downloadLink((end-4):end)==".sofa"
    
                targetfile = [targetDir '\' char(currentFiles(jj))];

                %% download file
                attempts = 0; % count download attempts
                for retry=1:downloadAttempts
                    try % to download the file
                        errorCatch=false;
                        websave(targetfile,downloadLink,options);
                    catch me
                        errorCatch=true;
                        if retry < downloadAttempts % give it another try
                            disp(['Download attempt failed, retry downloading (' num2str(retry) '/' num2str(downloadAttempts-1) ') from: ' downloadLink]);
                            appendToFile(logfile,[num2str(logIndex) char(9) num2str(ii) char(9) 'RETRY' char(9) 'Download attempt failed for the ' num2str(retry) '. time, wait and retry...' char(9) downloadLink]);   logIndex = logIndex + 1;
                            pause(retry * 5); % give it a break of a couple of seconds, server might block too many simulataneous downloads!?
                        else % ok, let's skip this file
                            warning(['Download' char(9) me.message]);
                        end
                    end
                    if errorCatch==false
                        break % leave loop
                    end
                end

                %% load file
                if errorCatch==false
                    try % to load the file
                        Obj=SOFAload(targetfile);
                    catch me
                        errorCatch=true;
                        warning(['SOFAload' char(9) me.message]);
                    end
                end

%                 % SOFAsave uses evalc instead of lastwarn
%                 [warnmsg, msgid] = lastwarn;
                if upgradeConventions == true
                    % upgrade SingleRoomDRIR2SingleRoomSRIR, MultiSpeakerBRIR2SingleRoomMIMOSRIR ?
                    switch Obj.GLOBAL_SOFAConventions
                        case 'SingleRoomDRIR'
                            Obj=SOFAconvertSingleRoomDRIR2SingleRoomSRIR(Obj);
                            warnmsg = 'Obj converted from SOFAconvertSingleRoomDRIR to SingleRoomSRIR convention.';  disp(warnmsg);
                            appendToFile(logfile,[num2str(logIndex) char(9) num2str(ii) char(9) 'WARNING' char(9) warnmsg char(9) targetfile]);   logIndex = logIndex + 1;
                        case 'MultiSpeakerBRIR'
                            Obj=SOFAconvertMultiSpeakerBRIR2SingleRoomMIMOSRIR(Obj);
                            warnmsg = 'Obj converted from MultiSpeakerBRIR to SingleRoomMIMOSRIR convention.';  disp(warnmsg);
                            appendToFile(logfile,[num2str(logIndex) char(9) num2str(ii) char(9) 'WARNING' char(9) warnmsg char(9) targetfile]);   logIndex = logIndex + 1;
                    end
                end
                %% save file
                if errorCatch==false
                    try % to save the file
%                         SOFAsave(temptargetfile,Obj,0); % save without compression (faster)
                        SOFAsaveResults = evalc('SOFAsave(temptargetfile,Obj,0);'); % evaluate function to get warnings
                        warnmsg = char(extractBetween(SOFAsaveResults,'Warning: ',']'));
%                         msgid = 'SOFA:save';
                        if ~isempty(warnmsg)
                            if size(warnmsg,1) > 1
                                warnmsg=replace(char(convertCharsToStrings(warnmsg')),'               ','; '); % multiple warning rows -> merge to one row
                            end
                            warning(['SOFAsave' char(9) warnmsg]);
                        end
                    catch me % error case
                        errorCatch=true;
                        warning(['SOFAsave' char(9) me.message]);
                    end
                end

                %% check for errors & warnings
                [warnmsg, msgid] = lastwarn;
                if ~isempty(warnmsg)
                    if errorCatch==true
                        appendToFile(logfile,   [num2str(logIndex)    char(9) num2str(ii) char(9) 'ERROR' char(9) warnmsg char(9) targetfile]);   logIndex = logIndex + 1;
                        appendToFile(logfileErr,[num2str(logErrIndex) char(9) num2str(ii) char(9) 'ERROR' char(9) warnmsg char(9) targetfile]);   logErrIndex = logErrIndex + 1;
                    else
                        appendToFile(logfile,[num2str(logIndex) char(9) num2str(ii) char(9) 'WARNING' char(9) warnmsg char(9) targetfile]);   logIndex = logIndex + 1;
                    end
                end
                lastwarn(''); errorCatch=false;
    
            end
            pause(0.01); % server errors occur after some files, if not paused here
        end
    else
        %% check for errors & warnings during webread & download
        [warnmsg, msgid] = lastwarn;
        if ~isempty(warnmsg)
            appendToFile(logfile,   [num2str(logIndex)     char(9) num2str(ii) char(9) 'ERROR' char(9) warnmsg char(9) currentSubLink]);   logIndex = logIndex + 1;
            appendToFile(logfileErr,[num2str(logErrIndex)  char(9) num2str(ii) char(9) 'ERROR' char(9) warnmsg char(9) currentSubLink]);   logErrIndex = logErrIndex + 1;
        end
        lastwarn(''); errorCatch=false;
    end
end


%% Epilogue
fid=fopen(logfile,'a+');     fprintf(fid, '\n\n%s', ['*** Checks done; end time: ' datestr(now(), 'dd.mm.yyyy - HH:MM:SS')]);
fid=fopen(logfileErr,'a+');  fprintf(fid, '\n\n%s', ['*** Checks done; end time: ' datestr(now(), 'dd.mm.yyyy - HH:MM:SS')]);

disp(['### See log file for errors and warnings: ' logfile])
disp(['### See log file for errors: ' logfileErr])
disp('      ');
disp('##############################################');
disp('##########   COMPLETED ALL CHECKS   ##########');
disp('##############################################');
toc; % timer

%% Helpers
function appendToFile(filename,text)
    fid = fopen(filename, 'a+');
    fprintf(fid, '\n%s', text);
    fclose(fid);
end
