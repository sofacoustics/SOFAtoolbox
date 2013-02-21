function [M,R,N,T] = SOFAcheckDimensions(Dataset)
%SOFACHECKDIMENSIONS 
%   [] = SOFAcheckDimensions(Dataset) checks dimensions of variables
%   contained in Dataset.

% SOFA API - function SOFAdisplay
% Copyright (C) 2012 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://www.osor.eu/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% -------------------------- prepare dimensions ---------------------------
VarNames = fieldnames(Dataset);
numVars = size(VarNames,1);
switch Dataset.DataType
	case 'FIR'
		[M R N] = size(Dataset.Data.FIR);
	case 'SpectralMagnitudePhase'
		[M R N] = size(Dataset.Data.Mag);
	otherwise
		error('Unknown DataType');
end
T = size(Dataset.TransmitterPosition,3);
SourceListenerVars=SOFAgetVariables('sourcelistener');
% TransmitterReceiverVars=SOFAgetVariables('transmitterreceiver');

% ------------------------ check matrix dimensions ------------------------
for ii=1:numVars
    CurrentValue = Dataset.(VarNames{ii});
    if strcmp(VarNames{ii},'Data') % do nothing (Data sets dimensions)
    elseif sum(strcmp(SourceListenerVars,VarNames{ii})) % Source/ListenerVars
        if ~((all(size(CurrentValue)==[M 3])) || (all(size(CurrentValue)==[1 3])))
            error('Dimensions of coordinate variables are not valid.');
        end
    elseif strcmp(VarNames{ii},'ReceiverPosition') % ReceiverPosition
        if ~((all(size(CurrentValue)==[1 3 1])) || (all(size(CurrentValue)==[M 3 1]) || ...
            (all(size(CurrentValue)==[1 3 R])) || (all(size(CurrentValue)==[M 3 R]))))
            error('Dimensions of ReceiverPosition are not valid.');
        end
    elseif strcmp(VarNames{ii},'TransmitterPosition') % TransmitterPosition
        % T doesn't need to be checked, as it is defined by the size of TransmitterPosition
        if ~((all(size(CurrentValue)==[1 3])) || (all(size(CurrentValue)==[M 3])))
            error('Dimensions of TransmitterPosition are not valid.');
        end
    elseif ~(size(size(CurrentValue),2)>2)
        % if matrix is not more than 2D -> check dimensions: [1 1], [M 1], [1 x], [M x]
        if ~(all(size(CurrentValue)==[1 1]) || all(size(CurrentValue)==[M 1]) || ...
            (size(CurrentValue,1)==1 && size(CurrentValue,2)>1) || ...
            (size(CurrentValue,1)==M && size(CurrentValue,2)>1))
            error('Invalid matrix dimensions.');
        end
    elseif size(size(CurrentValue),2)>2 % if matrix is more than 2D
        error('Invalid matrix dimensions.');
    end
end

end %of function