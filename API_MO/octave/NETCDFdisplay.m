function [] = NETCDFdisplay(filename)
%NETCDFDISPLAY 
%   [] = NETCDFdisplay(filename) displays information about specified SOFA file

% SOFA API - function octave/NETCDFdisplay
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

ncid=netcdf(filename,'r');

% ----- GLOBAL ATTRIBUTES ------------------------------------------------
globalAttr = ncatt(ncid);
fprintf('\n--- GLOBAL ATTRIBUTES ---\n');
for ii=1:length(globalAttr)
    attrName = ncname(globalAttr{ii});
    attrVal = globalAttr{ii}(:);
    fprintf('%s = %s\n',attrName,attrVal);
end


% ----- DIMENSIONS -------------------------------------------------------
dimensions = ncdim(ncid);
fprintf('\n--- DIMENSIONS ---\n');
for ii=1:length(dimensions)
    dimName = ncname(dimensions{ii});
    dimVal = dimensions{ii}(:);
    fprintf('%s = %s\n',dimName,num2str(dimVal));
end


% ----- VARIABLES --------------------------------------------------------
variables = ncvar(ncid);
fprintf('\n--- VARIABLES ---\n');
for ii=1:length(variables)
    varName = ncname(variables{ii});
    varVal = variables{ii}(:);
    fprintf('%s\n',varName);
    fprintf('\tSize: \t %s\n',num2str(size(varVal)));
    % get dimensions
    varDims = ncdim(variables{ii});
    varDimNames = [];
    for jj=1:length(varDims)
        varDimNames = [varDimNames ncname(varDims{jj})];
    end
    fprintf('\tDimensions: \t %s\n',varDimNames);
    % get attributes
    varAttr = ncatt(variables{ii});
    fprintf('\tAttributes:\n');
    for jj=1:length(varAttr)
        varAttrName = ncname(varAttr{jj});
        varAttrVal = varAttr{jj}(:);
        fprintf('\t\t%s = %s\n',varAttrName,varAttrVal);
    end
end

end %of function
