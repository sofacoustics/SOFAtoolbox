%COORDINATES Class useful to handl HRTF coordinates with different conventions
%   Usage: coordinates = Coordinates(SOFAobj);
%
%   Interface:
% 
%       Coordinates: constructor. The coordinates are converted into the 
%                    cartestian system is used.
%
%       return_positions: return the coordinates given a specific convention.
%
%       convert_positions: convert the stored coordinates given
%                    a specific convention.
%
%       find_position: return a position given the index of the coordinates
%                    matrix
%
%       normalize_distance: normalize the distance between receiver and
%                    source
% 
% 
%   Purpose:
%   Manage easily the coordinates system of the SOFA object.
% 
%   Requirements: 
%   -------------
%
%   1) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. thirdparty/SOFA)
% 
%
%   Examples:
%   ---------
% 
%   coordinates = Coordinates(SOFAobj);
%   horpolar = coords.return_positions('horizontal-polar');
%   
%   Load coordinates from SOFA object and convert them into the
%   horizontal-polar system. The output is organized as in the SOFA object.
%

% #Author: Roberto Barumerli
% #Author: Michael Mihocic: header documentation updated (28.10.2021)

classdef Coordinates < handle
    properties (SetAccess = private)
        pos(:,3) double {mustBeReal, mustBeFinite}
        pos_type string {mustBeMember(pos_type,{'horizontal-polar','spherical','cartesian', 'geodesic'})} = 'cartesian'
    end
    methods
        % constructor
        function obj = Coordinates(data, convention)
            if nargin == 1
                if strcmp(data.GLOBAL_Conventions, 'SOFA')
                    obj.pos = SOFAcalculateAPV(data);
                    obj.pos_type = data.SourcePosition_Type;
                else
                    error('Coordinates: SOFA object not valid!');
                end
            elseif nargin == 2
                if ismatrix(data)
                    obj.pos = data;
                    obj.pos_type = convention;
                else
                    error('Coordinates: position matrix not valid!');
                end
            else
                error('Coordinates: parameters not valid!');
            end
            
            obj.convert_positions('cartesian');
        end
        
        % return coordinates with a specific 
        function r = return_positions(obj, pos_type)
            if strcmp(pos_type, obj.pos_type)
                r = obj.pos;
            else
                r = SOFAconvertCoordinates(obj.pos,obj.pos_type, pos_type);
            end
        end
        
        function convert_positions(obj, pos_type)
            obj.pos = SOFAconvertCoordinates(obj.pos,obj.pos_type, pos_type);
            obj.pos_type = pos_type;
        end
        
        function r = count_positions(obj)
            r = size(obj.pos, 1);
        end
        
        function r = find_position(obj, idx)
            r = SOFAconvertCoordinates(obj.pos(idx,:),obj.pos_type, 'spherical');
        end
        
        function normalize_distance(obj)
            pos_type_temp = obj.pos_type;
            obj.convert_positions('spherical');
            obj.pos(:,3) = 1;
            obj.convert_positions(pos_type_temp);
        end
        
        function plot(obj)
            r = return_positions(obj, 'cartesian');
            figure
            scatter3(r(:,1),r(:,2),r(:,3),20,0.5*ones(size(r, 1),1),'filled');
            view([1 0 0])
            axis equal;   
        end

    end
    
	methods (Access = private)
    end
end