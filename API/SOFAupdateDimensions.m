function Obj = SOFAupdateDimensions(Obj)
%SOFAupdateDimensions
%   Obj = SOFAupdateDimensions(Obj) updates the dimensions in the SOFA
%   structure
%
%   Obj is a struct containing the data and meta.
%		The dimension variables are created and updated corresponding to the
%		conventions


% SOFA API - function SOFAupdateDimensions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or ñ as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
% You may not use this work except in compliance with the Licence.
% You may obtain a copy of the Licence at: http://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the Licence for the specific language governing  permissions and limitations under the Licence. 

%% General definitions
dims={'i';'r';'e';'n';'m';'c';'q'}; % dimensions
Obj.I=1;
OC = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'a');

%% Update dimension variables
f=fieldnames(rmfield(OC.Dimensions,'Data'));
for ii=2:length(dims)
	for jj=1:length(f)
		dim=strfind(OC.Dimensions.(f{jj}),dims{ii});
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.(upper(dims{ii}))=size(Obj.(f{jj}),dim(1));
% 			disp([dims{ii} ': ' f{jj} '(' num2str(dim(1)) ')']);
			break;
		end
	end
end

fd=fieldnames(OC.Dimensions.Data);
for ii=2:length(dims)
	for jj=1:length(fd)
		dim=strfind(OC.Dimensions.Data.(fd{jj}),dims{ii});
		if iscell(dim), dim=cell2mat(dim); end;
		if ~isempty(dim)
			Obj.(upper(dims{ii}))=size(Obj.Data.(fd{jj}),dim(1));
% 			disp([dims{ii} ': Data.' fd{jj} '(' num2str(dim(1)) ')']);
			break;
		end
	end
end
