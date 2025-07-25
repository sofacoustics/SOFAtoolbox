function Obj = SOFAupdateDimensions(Obj,varargin)
%SOFAupdateDimensions - Update the dimensions in the SOFA object
%   Usage: Obj = SOFAupdateDimensions(Obj)
%          Obj = SOFAupdateDimensions(Obj, flags)
%
%   Obj = SOFAupdateDimensions(Obj) updates the dimensions in the SOFA
%   object Obj. This can be done to check if all data and metadata comply
%   with the convention. It is automatically done by SOFAsave and SOFAload.
%
%   The updated dimension sizes are stored in Obj.API.X.
%
%   Obj = SOFAupdateDimensions(Obj, 'verbose',1) provides more information
%   on the update process, which can be usefull when Obj does not comply.
%
%   Obj = SOFAupdateDimensions(Obj, 'nodata') disables the checks on the
%   Data variables, which can be usefull when debugging the remaining
%   variables of Obj.

% #Author: Piotr Majdak
% #Author: Piotr Majdak: String support added (09.08.2014)
% #Author: Piotr Majdak: Verbose mode added (10.10.2020)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)
% #Author: Piotr Majdak: checkdim works now on string variables with singleton dimensions
% #Author: Michael Mihocic: Octave support fixed (29.11.2024)
%
% SOFA Toolbox - function SOFAupdateDimensions
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License.

definput.keyvals.Index=[];
definput.keyvals.verbose=0;
definput.flags.type={'data','nodata'};
[flags,kv]=SOFAarghelper({'Index'},definput,varargin);
v=kv.verbose;

%% Get conventions with allowed dimensions
OC = SOFAgetConventions(Obj.GLOBAL_SOFAConventions,'a');
if v, disp(['SOFA Convention: ' Obj.GLOBAL_SOFAConventions]); end

%% Add dimensions if required
dims=fieldnames(SOFAdefinitions('dimensions'));
for ii=1:size(dims,1)
    if ~isfield(Obj.API,dims{ii}), Obj.API.(dims{ii})=0; end
end

%% Update dimension sizes from the variables having dominant dimensions sizes
% fix dimension sizes
Obj.API.I=1;
Obj.API.C=3;
% check all metadata variables for dominant dimension sizes
dims='renm';
f=fieldnames(rmfield(OC.API.Dimensions,'Data'));
for ii=1:length(dims)
    for jj=1:length(f)
        dim=strfind(OC.API.Dimensions.(f{jj}),dims(ii));
        if iscell(dim), dim=cell2mat(dim); end;
        if ~isempty(dim)
            Obj.API.(upper(dims(ii)))=size(Obj.(f{jj}),dim(1));
            if (v), disp([upper(dims(ii)) ' set to ' num2str(size(Obj.(f{jj}),dim(1))) ' that is the #' num2str(dim(1)) ' dimension of ' f{jj}]); end
            break;
        end
    end
end
% check all data variables
if flags.do_data
    fd=fieldnames(OC.API.Dimensions.Data);
    for ii=1:length(dims)
        for jj=1:length(fd)
            dim=strfind(OC.API.Dimensions.Data.(fd{jj}),dims(ii));
            if iscell(dim), dim=cell2mat(dim); end;
            if ~isempty(dim)
                Obj.API.(upper(dims(ii)))=size(Obj.Data.(fd{jj}),dim(1));
                if (v), disp([upper(dims(ii)) ' set to ' num2str(size(Obj.Data.(fd{jj}),dim(1))) ' that is the #' num2str(dim(1)) ' dimension of Data.' fd{jj}]); end
                break;
            end
        end
    end
end

%% Update the dimensions of metadata variables
Smax=0;
X=rmfield(Obj,{'Data','API'});
if isfield(X,'PRIVATE'), X=rmfield(X,'PRIVATE'); end
Xf=fieldnames(X);
for ii=1:length(Xf)
    if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
        if isfield(OC.API.Dimensions, Xf{ii}), % is a known variable
            dim=OC.API.Dimensions.(Xf{ii});
            if ~iscell(dim), dim={dim}; end;
            [dim,S]=checkdim(Obj,dim,sizecell(Obj.(Xf{ii})));
            if isempty(dim),
                error([Xf{ii} ': dimension could not be matched.']);
            else
                Obj.API.Dimensions.(Xf{ii})=dim;
                if v, disp([Xf{ii} ': convention variable, used dimension: ' dim]); end;
            end
        else % is a user-defined variable
            if ~isfield(Obj.API.Dimensions,Xf{ii}),
                error([Xf{ii} ' seems to be a user-defined variable without dimension provided in API.Dimensions.']);
            else
                dim=Obj.API.Dimensions.(Xf{ii});
                [dim,S]=checkdim(Obj,{dim},sizecell(Obj.(Xf{ii})));
                if isempty(dim),
                    error([Xf{ii} ': dimension does not match.']);
                else
                    if v, disp([Xf{ii} ': user-defined variable, used dimension: ' dim]); end;
                end
            end
        end
        Smax=max(Smax,S);
    end
end
%% Update the dimensions of data variables
if flags.do_data
    Xf=fieldnames(Obj.Data);
    for ii=1:length(Xf)
        if isempty(strfind(Xf{ii},'_')),	% is not an attribute...
            if isfield(OC.API.Dimensions.Data, Xf{ii}), 			% is a known variable
                dim=OC.API.Dimensions.Data.(Xf{ii});
                if ~iscell(dim), dim={dim}; end;
                [dim,S]=checkdim(Obj,dim,sizecell(Obj.Data.(Xf{ii})));
                if isempty(dim),
                    error(['Data.' Xf{ii} ': dimension could not be matched.']);
                else
                    Obj.API.Dimensions.Data.(Xf{ii})=dim;
                    if v, disp(['Data.' Xf{ii} ': convention variable, used dimension: ' dim]); end;
                end
                Smax=max(Smax,S);
            else
                if ~isfield(Obj.API.Dimensions.Data,Xf{ii}),
                    error([Xf{ii} ' seems to be a user-defined variable without a dimension.']);
                else
                    dim=Obj.API.Dimensions.Data.(Xf{ii});
                    [dim,S]=checkdim(Obj,{dim},sizecell(Obj.Data.(Xf{ii})));
                    if isempty(dim),
                        error(['Data.' Xf{ii} ': dimension does not match.']);
                    else
                        if v, disp(['Data.' Xf{ii} ': user-defined variable, used dimension: ' dim]); end;
                    end
                end
            end
        end
    end
end
%% Update the size of the longest string
if Smax>0,
    Obj.API.S=Smax;
    if v, disp(['S set to ' num2str(Obj.API.S)]); end
else
    if v, disp('S unused (set to 0)'); end
end

%% Return the size of x. If x is a cell, return the size of the strings in x.
function s=sizecell(x,dim)
if iscell(x)
    s=size(char(x));
    if size(x,1)~=s(1) s=[size(x) s(2)]; end % multidim cellarays: s = [celldim1, celldim2, ... , celldimN, stringdim]
else
    s=size(x);
end

%% Get the sizes of the dimension variables according the dimension variables in str
function vec=getdim(Obj,str)
vec=arrayfun(@(f)(Obj.(f)),upper(str));

%% dims is a cell array with allowed dimensions.
% S is the size of the string dimension. S=0 when S does not exist
% dimA is a vector with the actual dimensions.
% dim is a string with the matching dimension
function [dim,S]=checkdim(Obj,dims,dimA)
dim=[]; S=0;
for jj=1:length(dims)
    dimS=dims{jj};
    if length(dimS)==1, dimS=[dimS 'I']; end % 1D required, but Matlab is always 2D at least.
    dimR=getdim(Obj.API,dimS);
    if ~isempty(strfind(dimS,'S')) % due to Octave support, do not replace by: if contains(dimS,'S')
        Sidx=strfind(dimS,'S');
        S=max(S,dimA(Sidx));
        dimR(Sidx)=dimA(Sidx); % string dim are always correct
    end
    if length(dimA)==length(dimR) % the same size?
        if dimA==dimR, dim=upper(dims{jj}); break; end	% found!
    elseif length(dimA)<length(dimR)	% extend the size?
        if [dimA ones(1,length(dimR)-length(dimA))]==dimR, dim=upper(dims{jj}); break; end; % found!
    end
end
