function results = SOFAget(Filename,MId,varargin)
%SOFAGETDATA
%   results = SOFAget(Filename,MId,varargin) reads several sets of data and Metadata from
%   a SOFA file for given IDs
%
%   Filename specifies the SOFA file from which the data is read.
%   MId is a vector of IDs specifying the measurement(s)
%   of which the data is read.
%   ReturnType is optional and specifies whether the function returns the
%   lodaded values as a struct or as a cell array. Default value is 'struct'.
%
%   If ReturnType is 'struct', the function returns a struct which contains
%   one field called 'Data' for the data and additional fields for each
%   metadata value. The name of these fields are identical to the names of the metadata.
%   If ReturnType is 'cell', the function returns a cell array with
%   the following structure:
%   results{x}{y}{id}
%   x ... number of variable
%   y = 1: variable name; y = 2: value
%   id ... adressing the different measurements

% SOFA API - function SOFAget
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% --------------- check and prepare variables ------------------
if(isnumeric(Filename))
  error('Filename must be a string.');
end
MId = MId(:);
if(~isnumeric(MId) | size(MId,2)>1)
  error('MId must be a numeric scalar or vector.');
end

ReturnType = 'struct'; % set default value for ReturnType
if(size(varargin,2)==1)
  varargin = cellstr(varargin);
  ReturnType = varargin{1};
end
if(isnumeric(ReturnType))
  error('ReturnType must be a string.');
end

if(strcmp(ReturnType,'struct'))
  results = struct; % initialize struct variable
elseif(strcmp(ReturnType,'cell'))
  % no need to initialize
else
  error('ReturnType must be either ''struct'' or ''cell''.');
end

%% ------------------------- get data ---------------------------
ncid = netcdf.open([Filename '.sofa'],'NC_NOWRITE');
try
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid); % get number of variables in file

for id=1:size(MId)
  count = 0;
  for ii=0:nvars-1 % loop through all variables in file
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,ii);
    % go through all dimensions of current variable and save their lengths to vectors
    for n=1:size(dimids,2)
      [DimName DimLength(n)] = netcdf.inqDim(ncid,dimids(n));
    end
    % ------ get Metadata ------
    if(strncmp(varname,'Data.',5))
      % do nothing (is done by SOFAgetData())
    elseif(size(dimids,2)==1) % scalar
      if(strcmp(ReturnType,'struct'))
        results.(varname){id} = netcdf.getVar(ncid,ii,0,1);
      elseif(strcmp(ReturnType,'cell'))
        results{count+1}{1} = varname;
        results{count+1}{2}{id} = netcdf.getVar(ncid,ii,0,1);
      end
    elseif(size(dimids,2)==2) % 2-D matrix
      if(DimLength(1)==1) % [1 x]
        if(strcmp(ReturnType,'struct'))
          results.(varname){id} = netcdf.getVar(ncid,ii,[0 0],[1 DimLength(2)]);
        elseif(strcmp(ReturnType,'cell'))
          results{count+1}{1} = varname;
          results{count+1}{2}{id} = netcdf.getVar(ncid,ii,[0 0],[1 DimLength(2)]);
        end
      elseif(DimLength(1)>1) % [M x]
        if(strcmp(ReturnType,'struct'))
          results.(varname){id} = netcdf.getVar(ncid,ii,[MId(id)-1 0],[1 DimLength(2)]);
        elseif(strcmp(ReturnType,'cell'))
          results{count+1}{1} = varname;
          results{count+1}{2}{id} = netcdf.getVar(ncid,ii,[MId(id)-1 0],[1 DimLength(2)]);
        end
      end
    elseif(size(dimids,2)==3) % 3-D matrix
      if(DimLength(1)==1) % [1 x y]
        if(strcmp(ReturnType,'struct'))
          results.(varname){id} = netcdf.getVar(ncid,ii,[0 0 0],[1 DimLength(2) DimLength(3)]);
        elseif(strcmp(ReturnType,'cell'))
          results{count+1}{1} = varname;
          results{count+1}{2}{id} = netcdf.getVar(ncid,ii,[0 0 0],[1 DimLength(2) DimLength(3)]);
        end
      elseif(DimLength(1)>1) % [M x y]
         if(strcmp(ReturnType,'struct'))
          results.(varname){id} = netcdf.getVar(ncid,ii,[MId(id)-1 0 0],[1 DimLength(2) DimLength(3)]);
        elseif(strcmp(ReturnType,'cell'))
          results{count+1}{1} = varname;
          results{count+1}{2}{id} = netcdf.getVar(ncid,ii,[MId(id)-1 0 0],[1 DimLength(2) DimLength(3)]);
         end
      end
    end
    if(~strncmp(varname,'Data.',5)) count = count + 1; end
  end
% --- get data ---
Data = SOFAgetData(Filename,MId,ReturnType);
% --- merge data ---
if(strcmp(ReturnType,'struct'))
  results.Data = Data;
elseif(strcmp(ReturnType,'cell'))
  results{count+1} = Data; % count has already been incremented after last Metadata was added
end
end % of for loop
catch
  if(exist('ncid','var') && ~isempty(ncid)) netcdf.close(ncid); end
  error(['An error occured during reading the SOFA file: ' lasterr()]);
  % TODO lasterr() should not be used any more...
end
netcdf.close(ncid);

end % of function