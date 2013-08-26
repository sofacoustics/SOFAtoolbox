function [Obj,Dims] = NETCDFload(filename,flags)
%NETCDFLOAD
%   [Obj,Dims] = NETCDFload(filename,'all') reads the SOFA object Obj with all
%       data from a SOFA file.
%
%   Obj = NETCDFload(filename,'nodata') ignores the Data variables while
%       reading.
%
%   Obj = NETCDFload(filename,[START COUNT]) reads only COUNT number of
%       measurements beginning with the index START.
%
%   [Obj,Dims] = NETCDFload(...) returns the dimension variables found in
%       the file as a string.

% SOFA API - function octave/NETCDFload
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 


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
        Obj.API.(fieldName) = fieldVal;
        dims{ii} = fieldName;
        startp.(fieldName) = 1;
        endp.(fieldName) = fieldVal;
    end
    Dims=cell2mat(dims)';

    % Check the requested measurements
    if isnumeric(flags)
        if Obj.API.M<flags(2)
            error('Requested end index exceeds the measurement count');
        end
        startp.M = flags(1);
        endp.M = flags(1)+flags(2)-1;
    end

    % ----- VARIABLES + ATTRIBUTES ---------------------------------------
    variables = ncvar(ncid);
    for ii=1:length(variables)
        fieldName = ncname(variables{ii});
        dims = ncdim(variables{ii});
        dimNames = [];
        for jj=1:length(dims)
            dimName = ncname(dims{jj});
            dimNames = [dimNames dimName];
        end
        % --- get data
        % check if we have something like Data.IR as fieldName and split it at
        % "."
        if strfind(fieldName,'Data.')
            if ~strcmp(flags,'nodata')
                fieldName1 = fieldName(1:strfind(fieldName,'.')-1);
                fieldName2 = fieldName(strfind(fieldName,'.')+1:end);
                Obj.API.Dimensions.(fieldName1).(fieldName2) = dimNames;
                if length(dimNames)==3
                    Obj.(fieldName1).(fieldName2) = ...
                        variables{ii}(startp.(dimNames(1)):endp.(dimNames(1)),...
                                      startp.(dimNames(2)):endp.(dimNames(2)),...
                                      startp.(dimNames(3)):endp.(dimNames(3)));
                elseif length(dimNames)==2
                    Obj.(fieldName1).(fieldName2) = ...
                        variables{ii}(startp.(dimNames(1)):endp.(dimNames(1)),...
                                      startp.(dimNames(2)):endp.(dimNames(2)));
                else
                    Obj.(fieldName1).(fieldName2) = ...
                        variables{ii}(startp.(dimNames(1)):endp.(dimNames(1)));
                end
            end
        else
            Obj.API.Dimensions.(fieldName) = dimNames;
            if length(dimNames)==3
                Obj.(fieldName) = ...
                    variables{ii}(startp.(dimNames(1)):endp.(dimNames(1)),...
                                  startp.(dimNames(2)):endp.(dimNames(2)),...
                                  startp.(dimNames(3)):endp.(dimNames(3)));
            elseif length(dimNames)==2
                Obj.(fieldName) = ...
                    variables{ii}(startp.(dimNames(1)):endp.(dimNames(1)),...
                                  startp.(dimNames(2)):endp.(dimNames(2)));
            else
                Obj.(fieldName) = ...
                    variables{ii}(startp.(dimNames(1)):endp.(dimNames(1)));
            end
        end
           
        % --- get attributes
        attr = ncatt(variables{ii});
        for jj=1:length(attr)
            attrName = ncname(attr{jj});
            attrVal = attr{jj}(:);
            % check if we have something like Data.IR as fieldName and split it at
            % "."
            if strfind(fieldName,'Data.')
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
    error(['An error occured during reading the SOFA file: ' ...
        lasterror.message lasterror.stack]);
end_try_catch

close(ncid)

end %of function
