function [varName,varContent] = NETCDFload(Filename,varargin)
%NETCDFLOAD
%   results = NETCDFload(Filename,ReturnType) reads all data (no metadata) from
%   a SOFA file.
%
%   Filename specifies the SOFA file from which the data is read.

% SOFA API - function SOFAloadData
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences;
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

if isempty(varargin)
    varargin={'all'};
end

%% --------------------------- N E T C D F load ---------------------------
try
    ncid = netcdf.open(Filename,'NC_NOWRITE');% open file
    [~,numVars,~,~] = netcdf.inq(ncid); % get number of variables stored in file
    
    count=1;
    for ii=1:numVars % LOOP through all variables in file
        currentVarName = netcdf.inqVar(ncid,ii-1); % get current variable name
        for jj=1:length(varargin)
            switch varargin{jj}
                case 'all'
                    varName{count}=currentVarName;
                    varContent{count}=netcdf.getVar(ncid,ii-1); % get values from current variable
                    count=count+1;
                case 'meta'
                    if ~strncmp(currentVarName,'Data.',5) % if current variable is Metadata
                        varName{count}=currentVarName;
                        varContent{count}=netcdf.getVar(ncid,ii-1); % get values from current variable
                        count=count+1;
                    end
                case 'data'
                    if strncmp(currentVarName,'Data.',5) % if current variable is Metadata
                        varName{count}=currentVarName;
                        varContent{count}=netcdf.getVar(ncid,ii-1); % get values from current variable
                        count=count+1;
                    end
                case currentVarName
                    varName{count}=currentVarName;
                    varContent{count}=netcdf.getVar(ncid,ii-1); % get values from current variable
                    count=count+1;
                otherwise
                    error('Wrong SOFA variable name!')
            end
        end
    end
    netcdf.close(ncid)
catch netcdfError
    netcdf.abort(ncid);
    throw(netcdfError)
end

end %of function