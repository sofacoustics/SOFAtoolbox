%test_SOFAall - Test script, running all demos from 'demos' subfolder and testing some SOFA Toolbox functionalities.

% #Author: Michael Mihocic: header documentation updated (28.10.2021)
%
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% Preparations
clc; close all; clear; % clean-up first
tic; % timer
SOFAstart('restart');
% warning('off','SOFA:upgrade');
% warning('off','SOFA:load');
% warning('off','SOFA:save');

%% Prologue
disp('      ');
disp('############################################');
disp('##########   TEST SOFA URL DATA   ##########');
disp('############################################');
disp('      ');
% disp('!!! Make sure that all source files are available in the data subdirectories. See individual readme.txt files for more information !!!');
% disp('      ');

%% log file
logfile=[SOFAdbPath '\urlDatabase\log.txt'];
% need to create sub directory first?
if ~exist([SOFAdbPath '\urlDatabase\'], 'dir')
   mkdir([SOFAdbPath '\urlDatabase\']);
end
fid=fopen(logfile,'w');
fprintf(fid, '%s\n',[ '*** Checking SOFA files for errors and warnings while downloading, loading & saving; start time: ' char(datetime)]);
fclose(fid);
lastwarn(''); % clear last warning
errorCatch=false; % use variable to detect catched errors

%% Get subdirectories from URL db
dirDatabase=[SOFAdbURL '/database/'];
options = weboptions("timeout", 60);
dirContent=webread(dirDatabase,options);

% Get subdirectories in database
subDirs = extractBetween(dirContent, "alt=""[DIR]""></td><td><a href=""", "/""");

%% Download data
for ii=1:size(subDirs,1)
    errorCatch=false; % ok, give it a chance

    % let's scan the current subdirectory
    currentSubLink=[dirDatabase char(subDirs(ii))];
    options = weboptions("timeout", 60);
    try
        currentContent=webread(currentSubLink,options);
    catch me
        errorCatch=true;
        warning([me.message]);
        pause(5); % give it a break, server might block multiple downloads!?
    end

    if errorCatch==false
        % let's look for files in this subdirectory
        currentFiles = extractBetween(currentContent, "</td><td><a href=""", """>");
        % where do we wanna store the data?
        targetDir=[SOFAdbPath '\urlDatabase\' char(strrep(subDirs(ii),"%20", " "))];
    
        % a short status message:
        disp(['### Downloading URL data: ' char(strrep(subDirs(ii),"%20", " "))]);
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
                try % to download the file
                    websave(targetfile,downloadLink,options);
                catch me
                    errorCatch=true;
                    warning(['Download file: ' me.message]);
                    pause(5); % give it a break, server might block multiple downloads!?
                end

                %% load file
                if errorCatch==false
                    try % to load the file
                        Obj=SOFAload(targetfile);
                    catch me
                        errorCatch=true;
                        warning(['Load file: ' me.message]);
                    end
                end

                %% save file
                if errorCatch==false
                    try % to save the file
                        targetfile = [targetDir '\' strrep(char(currentFiles(jj)),'.sofa','_loadsave.sofa')];
                        SOFAsave(targetfile,Obj);
                    catch me
                        errorCatch=true;
                        warning(['Save file: ' me.message]);
                    end
                end

                %% check for errors & warnings
                [warnmsg, msgid] = lastwarn;
                if ~isempty(warnmsg)
                    appendToFile(logfile,[targetfile ': ' warnmsg])
                end
                lastwarn(''); errorCatch=false;
    
            end
            pause(0.1); % server errors occur after some files, if not paused here
        end
    else
        %% check for errors & warnings during download
        [warnmsg, msgid] = lastwarn;
        if ~isempty(warnmsg)
            appendToFile(logfile,[currentSubLink ': ' warnmsg])
        end
        lastwarn(''); errorCatch=false;
    end
end



%% Epilogue

fid=fopen(logfile,'a+');
fprintf(fid, '\n\n%s', ['*** Checks done; end time: ' char(datetime)]);

disp('##############################################');
disp('####   COMPLETED ALL CHECKS SUCCESSFULLY  ####');
disp('##############################################');
toc; % timer

function appendToFile(filename,text)
    fid = fopen(filename, 'a+');
    fprintf(fid, '\n%s', text);
    fclose(fid);
end
