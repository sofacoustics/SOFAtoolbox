function [varName,varContent] = NETCDFload(filename,varargin)
%NETCDFLOAD
%   [varName,varContent] = NETCDFload(filename,ReturnType) reads all data (no metadata) from
%   a SOFA file.
%
%   filename specifies the SOFA file from which the data is read.

% SOFA API - function octave/NETCDFload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

if isempty(varargin)
    varargin={'all'};
end

%% --------------------------- N E T C D F load ---------------------------
try
    ncid = netcdf(filename,'r','netcdf4'); % open file
    sofa = ncvar(ncid); % get the variables stored in the file
    numVars = length(sofa); % get number of variables in file    

    count=1;
    loadVar=0;
    for ii=1:numVars % LOOP through all variables in file
        currentVarName = ncname(sofa{ii}); % get current variable name
        for jj=1:length(varargin)
            switch varargin{jj}
                case 'all'
                    loadVar=1;
                case 'meta'
                    if ~strncmp(currentVarName,'Data.',5) % if current variable is Metadata
                        loadVar=1;
                    end
                case 'data'
                    if strncmp(currentVarName,'Data.',5) % if current variable is Data
                        loadVar=1;
                    end
                case currentVarName
                    loadVar=1;
                otherwise
                    error('Wrong SOFA variable name!')
            end
            if loadVar
                varName{count}=currentVarName;
                currentVarValue=ncid{currentVarName}(:);
                if length(ncdim(sofa{ii})(:))<3
                    varContent{count}=transpose(currentVarValue);
                elseif length(ncdim(sofa{ii})(:))==3
                    varContent{count}=permute(currentVarValue,[3 2 1]);
                end
                count=count+1;
                loadVar=0;
            end
        end
    end
catch
    if exist('ncid','var') && ~isempty(ncid)
        netcdf.close(ncid);
    end
    error(['An error occured during reading the SOFA file:\n' lasterror.message]);
end_try_catch

close(ncid)

end %of function