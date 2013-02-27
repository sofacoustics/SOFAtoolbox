function DimNames = SOFAgetDimensions()
%SOFAGETDIMENSIONS
%   DimNames = SOFAgetDimensions() returns a struct with all dimension names

% SOFA API - function SOFAgetDimensions
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% ----- Data matrix ------------------------------------------------------
DimNames.Measurements = 'Measurements';
DimNames.Receivers = 'Receivers';
DimNames.Samples = 'Samples';
%DimNames.Transmitters = 'Transmitters';

%
DimNames.Coordinates = 'Coordinates';
DimNames.Scalar = 'Scalar';
DimNames.Unlimited = 'Unlimited';
%DimNames.String = 'String';
%DimNames.Vector = 'Vector';

end %of function