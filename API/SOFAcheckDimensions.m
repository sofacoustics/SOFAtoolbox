function [M,R,N,T] = SOFAcheckDimensions(dataset)
%SOFACHECKDIMENSIONS 
%   [] = SOFAcheckDimensions(dataset) checks dimensions of variables
%   contained in dataset.

% SOFA API - function SOFAdisplay
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

% -------------------------- prepare dimensions ---------------------------
varNames = fieldnames(dataset);
numVars = size(varNames,1);
switch dataset.DataType
	case 'IR'
		[M, R, N] = size(dataset.Data.IR);
	case 'SpectralMagnitudePhase'
		[M, R, N] = size(dataset.Data.Mag);
	otherwise
		error('Unknown DataType');
end
T = size(dataset.TransmitterPosition,1);
% SourceListenerVars=SOFAgetVariables('sourcelistener');
% TransmitterReceiverVars=SOFAgetVariables('transmitterreceiver');

% ------------------------ check matrix dimensions ------------------------
for ii=1:numVars
    currentValue = dataset.(varNames{ii});
    if strcmp(varNames{ii},'Data') % do nothing (Data sets dimensions)
    elseif strcmp(currentValue,'')
%     elseif sum(strcmp(SourceListenerVars,varNames{ii})) % Source/ListenerVars
%         if ~((all(size(currentValue)==[M 3])) || (all(size(currentValue)==[1 3])))
%             error(['Dimensions of coordinate variables are not valid: ' num2str(ii)]);
%         end
    elseif strcmp(varNames{ii},'ReceiverPosition') % ReceiverPosition
        if ~((all(size(currentValue)==[R 3])) || (all(size(currentValue)==[R 3 M])))
            error('Dimensions of ReceiverPosition are not valid.');
        end
    elseif strcmp(varNames{ii},'TransmitterPosition') % TransmitterPosition
        % T doesn't need to be checked, as it is defined by the size of TransmitterPosition
        if ~((all(size(currentValue)==[T 3])) || (all(size(currentValue)==[T M 3])))
            error('Dimensions of TransmitterPosition are not valid.');
        end
    elseif ~(size(size(currentValue),2)>2)
        % if matrix is not more than 2D -> check dimensions: [1 1], [M 1], [1 x], [M x]
        if ~(all(size(currentValue)==[1 1]) || all(size(currentValue)==[M 1]) || ...
            (size(currentValue,1)==1 && size(currentValue,2)>1) || ...
            (size(currentValue,1)==M && size(currentValue,2)>1))
            error('Invalid matrix dimensions.');
        end
    elseif size(size(currentValue),2)>2 % if matrix is more than 2D
        error('Invalid matrix dimensions.');
    end
end

end %of function