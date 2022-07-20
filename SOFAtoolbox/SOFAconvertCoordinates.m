 function output = SOFAconvertCoordinates(input,input_type,output_type,~,~)
%SOFAconvertCoordinates
%   output = SOFAconvertCoordinates(input,input_type,output_type), converts the specified coordinate variable to specified
%   output_type and returns the results.
% 

% #Author: Piotr Majdak
% #Author: Michael Mihocic: type horizontal-polar removed (not defined in SOFA) (08.03.2021)
% #Author: Michael Mihocic: doc fixed, header documentation updated (28.10.2021)
% 
% SOFA API - function SOFAconvertCoordinates
% Copyright (C) 2012-2021 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

%% check input
if strcmp(input_type,'cartesian')==0 && ...
        strcmp(input_type,'spherical')==0 && ...
        strcmp(input_type,'geodesic')==0 % && ...
%         strcmp(input_type,'horizontal-polar')==0
    error('Specified "input_type" is not supported');
end
if strcmp(output_type,'cartesian')==0 && ...
        strcmp(output_type,'spherical')==0 && ...
        strcmp(output_type,'geodesic')==0 % && ...
%         strcmp(output_type,'horizontal-polar')==0
    error('Specified "output_type" is not supported');
end

output=input;
%% convert coordinates if necessary
if strcmp(output_type,input_type)==0
    temp=input;
    switch input_type
        case 'cartesian'
            %do nothing
        case {'spherical','geodesic'}
            [temp(:,1),temp(:,2),temp(:,3)]=sph2cart(deg2rad(input(:,1)),deg2rad(input(:,2)),input(:,3));
%         case 'horizontal-polar'
%             [temp(:,1),temp(:,3),temp(:,2)]=sph2cart(deg2rad(input(:,2)),deg2rad(input(:,1)),input(:,3));
%             temp(:,2)=-temp(:,2);
    end

    output=temp;
    switch output_type
        case 'cartesian'
            %do nothing
        case {'spherical','geodesic'}
            [output(:,1),output(:,2),output(:,3)]=cart2sph(temp(:,1),temp(:,2),temp(:,3));
            output(:,1:2)=rad2deg(output(:,1:2));
%         case 'horizontal-polar'
% %             [output(:,2),output(:,1),output(:,3)]=cart2sph(temp(:,1),temp(:,3),-temp(:,2));
% %             output(:,1:2)=rad2deg(output(:,1:2));
% %             output(:,1)=-output(:,1);
% %             output(:,2)= mod(output(:,2)+90,360)-90;
%         [output(:,1),output(:,2),output(:,3)]=cart2sph(temp(:,1),temp(:,3),-temp(:,2));
%         [output(:,1),output(:,2)]=sph2hor(output(:,1),output(:,2));
    end
end


    
    