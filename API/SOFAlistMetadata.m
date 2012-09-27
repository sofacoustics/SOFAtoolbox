function [] = SOFAlistMetadata(Filename)
%SOFALISTVARIABLES 
%   [] = SOFAlistMetadata(Filename) lists all Metadata contained in a SOFA file.
% 
%	The function opens the SOFA file, and prints the names of all Metadata 
%   variables contained in the SOFA file specified by Filename.

% SOFA API - function SOFAlistVariables
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences; Wolfgang Hrauda
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

if(isnumeric(Filename))
  error('Filename must be a string.');
end

Metadata = SOFAloadMetadata(Filename,'struct');
fprintf(['Metadata in file ' Filename '.sofa:\n\n']);
disp(fieldnames(Metadata));

end % of function