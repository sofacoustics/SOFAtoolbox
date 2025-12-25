function [output,output_units] = SOFAconvertCoordinates(input,input_type,output_type,input_units,~)
%SOFAconvertCoordinates - Convert between the coordinates systems
%   Usage: output = SOFAconvertCoordinates(input,input_type,output_type)
%     [output,output_units] = SOFAconvertCoordinates(input,input_type,output_type,input_units)
%
%   SOFAconvertCoordinates(input,input_type,output_type), converts input given in 
%   the coordinate system input_type to the coordinate system described by output_type.
%   input_type and output_type can be 'cartesian' or 'spherical' as specified in AES69.
%   input must be a matrix of X-by-C.
%
%   [output,output_units] = SOFAconvertCoordinates(...,input_units)
%   converts the units as well. 

% #Author: Piotr Majdak
% #Author: Michael Mihocic: type horizontal-polar removed (not defined in SOFA) (08.03.2021)
% #Author: Michael Mihocic: doc fixed, header documentation updated (28.10.2021)
% #Author: Piotr Majdak: convertion of units added (25.12.2025)
% 
% SOFA Toolbox - function SOFAconvertCoordinates
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% check the input type
input_type=lower(input_type);
if strcmp(input_type,'cartesian')==0 && ...
        strcmp(input_type,'spherical')==0 && ...
        strcmp(input_type,'geodesic')==0 
    error('Specified "input_type" is not supported');
end
output_type=lower(output_type);
if strcmp(output_type,'cartesian')==0 && ...
        strcmp(output_type,'spherical')==0 && ...
        strcmp(output_type,'geodesic')==0 
    error('Specified "output_type" is not supported');
end

output=input;
if ~exist('input_units','var')
  input_units='unknown,unknown,unknown'; 
  output_units = input_units;
else
    % fix units if required
  parts = strsplit(input_units, ',');
  if length(parts)==1
    output_units = [parts{1} ',' parts{1} ',' parts{1}]; % only a single dimension provided, use that for all dimensions
  else
    output_units = [parts{1} ',' parts{2} ',' parts{3}]; % preserve units for all three dimensions
  end
end

%% convert to Cartesian (if necessary) and then to the output type
if strcmp(output_type,input_type)==0
    temp=input;
    switch input_type
        case 'cartesian'
           

        case {'spherical','geodesic'}
            [temp(:,1),temp(:,2),temp(:,3)]=sph2cart(deg2rad(input(:,1)),deg2rad(input(:,2)),input(:,3));
            parts = strsplit(input_units, ',');
            output_units = [parts{3} ',' parts{3} ',' parts{3}]; % use the unit of the radius
    end

    output=temp;
    switch output_type
        case 'cartesian'
            %do nothing
        case {'spherical','geodesic'}
            [output(:,1),output(:,2),output(:,3)]=cart2sph(temp(:,1),temp(:,2),temp(:,3));
            output(:,1:2)=rad2deg(output(:,1:2));
            parts = strsplit(output_units, ',');
            output_units = ['degree,degree,' parts{1}]; 
    end
end


    
    