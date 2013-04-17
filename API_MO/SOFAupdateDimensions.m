function [Obj,dims] = SOFAupdateDimensions(Obj)
%SOFAupdateDimensions
%   [Obj,dims] = SOFAupdateDimensions(Obj) updates the dimensions in the SOFA
%   structure
%
%   Obj is a struct containing the data and meta.
%		The dimension variables are created and updated corresponding to the
%		conventions
%   dims is a struct containing all dimension variables

% SOFA API - function SOFAupdateDimensions
% Copyright (C) 2012-2013 Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.1 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence")
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
			break;
		end
	end
end

%% create field names which should have dimensions
X=rmfield(rmfield(Obj,'Data'),'Dimensions');
for ii=1:length(dims)
	if isfield(X,upper(dims{ii})), X=rmfield(X,upper(dims{ii})); end
end
Xf=fieldnames(X);

%% Update the dimensions structure w/o data
for ii=1:length(Xf)
	if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
		if isfield(OC.Dimensions, Xf{ii}), % is a known variable		
			dim=OC.Dimensions.(Xf{ii});
			if ~iscell(dim), dim={dim}; end;
			dim=checkdim(Obj,dim,size(Obj.(Xf{ii})));
			if isempty(dim),
				error([Xf{ii} ': dimension could not be matched.']);
			else
				Obj.Dimensions.(Xf{ii})=dim;
			end
		else % is a user-defined variable						
			if ~isfield(Obj,Xf{ii}),
				error([Xf{ii} ' seems to be a user-defined variable without a dimension.']);
			end
		end		
	end
end

%% Update/check dimensions of data
Xf=fieldnames(Obj.Data);
for ii=1:length(Xf)
	if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
		if isfield(OC.Dimensions.Data, Xf{ii}), 			% is a known variable
			dim=OC.Dimensions.Data.(Xf{ii}); 
			if ~iscell(dim), dim={dim}; end;
			dim=checkdim(Obj,dim,size(Obj.Data.(Xf{ii})));
			if isempty(dim),
				error(['Data.' Xf{ii} ': dimension could not be matched.']);
			else
				Obj.Dimensions.Data.(Xf{ii})=dim;
			end
		else
			error(['Unknown data variable' Xf{ii} '.']);
		end		
	end
end

%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
	vec=arrayfun(@(f)(Obj.(f)),upper(str));

%% dims is a cell array with allowed dimensions. dimA is a vector with the actual dimensions.
% dim is a string with the matching dimension
function dim=checkdim(Obj,dims,dimA)
dim=[];
if prod(dimA)==1, dimA=1; end; % squeeze scalars
for jj=1:length(dims)
	dimR=getdim(Obj,dims{jj});
	if length(dimA)==length(dimR), % the same size?
		if dimA==dimR, dim=upper(dims{jj}); break; end;	% found!
	elseif length(dimA)<length(dimR)	% extend the size?
		if [dimA ones(1,length(dimR)-length(dimA))]==dimR, dim=upper(dims{jj}); break; end; % found!
	end
end