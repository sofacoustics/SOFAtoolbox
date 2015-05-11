function SOFAcompileConventions(conventions)
%SOFAcompileConventions
%
%   Obj = SOFAcompileConventions(sofaconventions) compiles the specified
%   SOFA conventions. For every convention a CSV file has to exist which
%   will be compiled to a .mat file used later by SOFAgetConventions().
% 
%   The CSV file must be in the directory conventions and have the same
%   filename as conventions. SOFAcompileConventions generates 3 files, one
%   for each flag (r, m, and all). 
%
%   Before compiling, SOFAcompileConventions checks if the modification
%   date of the .mat files is older than that of the .csv file. Compiling
%   is not performed if all .mat files are newer than the .csv file. This
%   behaviour is required for operation in a read-only directory. 

% SOFA API 
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

baseFolder = fileparts(which('SOFAstart'));

if nargin<1
    conventionFiles = dir(fullfile(baseFolder,'conventions','*.csv'));
    conventions={};
    for file = conventionFiles'
        [~,name,ext] = fileparts(file.name);
        % Check if mat files exist for every convention flag (r,m,a)
        flagsCounter = 0;
        for flag = 'rma'
            flagFile = dir(fullfile(baseFolder,'conventions', ...
                             strcat(name,'-',flag,'.mat')));
            if ~isempty(flagFile) && flagFile(1).datenum>file.datenum
                flagsCounter = flagsCounter+1;
            end
        end
        % If not all three files are up to request conventions compilation
        if flagsCounter~=3
            conventions{end+1} = name;
        end
    end
elseif ~iscell(conventions)
    conventions={conventions};
end

  
%% ----- Convert convention csv files into mat files -----
for convention = conventions
    % Read convention description from csv file
    fid = fopen(fullfile(baseFolder,'conventions', ...
                         strcat(convention{:},'.csv')));
    C = textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','Headerlines',1);
    fclose(fid);
    % Convert to mat files for r,m,a cases
    for flag = 'rma'
        % Convert to SOFA object
        Obj = compileConvention(C,flag);
        % Write to mat file
        if strcmp(Obj.GLOBAL_SOFAConventions,convention{:})
            if strcmp(flag,'r') % Display message only the very first time
                disp(['Compiling ',convention{:},' ', ...
                            Obj.GLOBAL_SOFAConventionsVersion]);
            end
            save(fullfile(baseFolder,'conventions', ...
                          strcat(convention{:},'-',flag,'.mat')), ...
                 'Obj','-v7');
        end
    end
end
end % of main function


%% ----- Subroutines -----------------------------------------------------
function Obj = compileConvention(convention,flag)
    % Compile convention mat structure for the specified flag
    %
    % The csv files provide the following columns (corresponding cell numbers in
    % brackets)
    % Name {1}, Default {2}, Flags {3}, Dimensions {4}, Type {5}, Comment {6}
    convName = convention{1};
    convDefault = convention{2};
    convFlags = convention{3};
    convDimensions = convention{4};
    convType = convention{5};
    convComment = convention{6};

    % Create object structure
    for ii=1:length(convName)
        % Append 'a' to Flags entry as it only contains 'm' or 'r' in the csv file
        convFlags{ii} = strcat(convFlags{ii},'a');
        if ~isempty(regexp(convFlags{ii},flag))
            var = regexprep(convName{ii},':','_');
            switch lower(convType{ii})
            case 'double'
                % Convert default to double
                convDefault{ii} = str2num(convDefault{ii});
            case 'string'
                eval(['convDefault{ii}=' convDefault{ii} ';']);
            end
            if isempty(strfind(var,'Data.'))
                Obj.(var) = convDefault{ii};
                if isempty(strfind(var,'_')) % && ~sum(strcmp(var,dims))
                    x2 = regexprep(convDimensions{ii},' ',''); %  remove spaces
                    y = regexprep(x2,',',['''' 10 '''']); % enclose in quotations and insert line breaks
                    Obj.API.Dimensions.(var)=eval(['{''' y '''}']);
                end
            else      
                Obj.Data.(var(6:end)) = convDefault{ii};
                if isempty(strfind(var(6:end),'_')) 
                    x2 = regexprep(convDimensions{ii},' ',''); %  remove spaces
                    y = regexprep(x2,',',['''' 10 '''']); % enclose in quatations and insert line breaks
                    Obj.API.Dimensions.Data.(var(6:end))=eval(['{''' y '''}']);
                end      
            end
        end
    end


    % ----- Overwrite some special fields -----
    if isfield(Obj,'GLOBAL_APIVersion')
        Obj.GLOBAL_APIVersion = SOFAgetVersion;
    end
    if isfield(Obj,'GLOBAL_APIName')
        Obj.GLOBAL_APIName = 'ARI Matlab/Octave API';
    end

    % ----- Create dimension size variables - if not read-only -----
    if strcmp(flag,'r')
        return;
    else
        % Fix dimension sizes (why we have to fix them?)
        Obj.API.I = 1;
        Obj.API.C = 3;
        % Variable-dependent dimension sizes
        dims = 'renm';
        % Check all metadata variables
        fields =fieldnames(rmfield(Obj.API.Dimensions,'Data'));
        for ii=1:length(dims)
            for jj=1:length(fields)
                dim = strfind(Obj.API.Dimensions.(fields{jj}),dims(ii));
                if iscell(dim), dim=cell2mat(dim); end;
                if ~isempty(dim)
                    Obj.API.(upper(dims(ii)))=size(Obj.(fields{jj}),dim(1));
                    break;
                end
            end
        end
        % Check all data variables
        fields = fieldnames(Obj.API.Dimensions.Data);
        for ii=1:length(dims)
            for jj=1:length(fields)
                dim = strfind(Obj.API.Dimensions.Data.(fields{jj}),dims(ii));
                if iscell(dim), dim=cell2mat(dim); end;
                if ~isempty(dim)
                    Obj.API.(upper(dims(ii)))=size(Obj.Data.(fields{jj}),dim(1));
                    break;
                end
            end
        end
    end
end
