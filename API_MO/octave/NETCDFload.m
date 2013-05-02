function Obj = NETCDFload(filename,flags)
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


glob='GLOBAL_';

%% --------------------------- N E T C D F load --------------------------
try
    ncid = netcdf(filename,'r','netcdf4'); % open file
    Obj = [];
    
    % ----- GLOBAL ATTRIBUTES --------------------------------------------
    globalAttr = ncatt(ncid);
    for ii=1:length(globalAttr)
        fieldName = ncname(globalAttr{ii});
        fieldVal = globalAttr{ii}(:);
        Obj.([glob fieldName]) = fieldVal;
    end

    % ----- DIMENSIONS ---------------------------------------------------
    dimensions = ncdim(ncid);
    for ii=1:length(dimensions)
        fieldName = ncname(dimensions{ii});
        fieldVal = dimensions{ii}(:);
        Obj.(fieldName) = fieldVal;
        dims{ii} = fieldName;
        startp(ii) = 0;
        countp(ii) = fieldVal;
    end
    Dims=cell2mat(dims)';

    %% Check the requested measurements
    %if isnumeric(flags)
    %    if Obj.M<flags(2)
    %        error('Requested end index exceeds the measurement count');
    %    end
    %    startp(strfind(Dims,'M'))=flags(1)-1;
    %    countp(strfind(Dims,'M'))=flags(2);
    %end

    % ----- VARIABLES + ATTRIBUTES ---------------------------------------
    variables = ncvar(ncid);
    for ii=1:length(variables)
        fieldName = ncname(variables{ii});
        fieldVal = variables{ii}(:);
        if isempty(strfind(Dims,fieldName)) % don't store data and dimensions for dimensions
            % --- get data
            % check if we have something like Data.IR as fieldName and split it at
            % "."
            if strfind(fieldName,'.')
                fieldName1 = fieldName(1:strfind(fieldName,'.')-1);
                fieldName2 = fieldName(strfind(fieldName,'.')+1:end);
                Obj.(fieldName1).(fieldName2) = fieldVal;
            else
                Obj.(fieldName) = fieldVal;
            end
        
            % --- get dimensions
            dims = ncdim(variables{ii});
            dimNames = [];
            for jj=1:length(dims)
                dimName = ncname(dims{jj});
                dimNames = [dimNames dimName];
            end
            % check if we have something like Data.IR as fieldName and split it at
            % "."
            if strfind(fieldName,'.')
                fieldName1 = fieldName(1:strfind(fieldName,'.')-1);
                fieldName2 = fieldName(strfind(fieldName,'.')+1:end);
                Obj.Dimensions.(fieldName1).(fieldName2) = dimNames;
            else
                Obj.Dimensions.(fieldName) = dimNames;
            end
        end
            
        % --- get attributes
        attr = ncatt(variables{ii});
        for jj=1:length(attr)
            attrName = ncname(attr{jj});
            attrVal = attr{jj}(:);
            % check if we have something like Data.IR as fieldName and split it at
            % "."
            if strfind(fieldName,'.')
                fieldName1 = fieldName(1:strfind(fieldName,'.')-1);
                fieldName2 = fieldName(strfind(fieldName,'.')+1:end);
                Obj.(fieldName1).([fieldName2 '_' attrName]) = attrVal;
            else
                Obj.([fieldName '_' attrName]) = attrVal;
            end
        end
          
    end

catch
    if exist('ncid','var') && ~isempty(ncid)
        close(ncid);
    end
    error(['An error occured during reading the SOFA file: ' lasterror.message]);
end_try_catch

close(ncid)

end %of function
