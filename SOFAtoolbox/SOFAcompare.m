function [tf,reason,where] = SOFAcompare(Obj1, Obj2, varargin)
%SOFAcompare - compare two SOFA objects
%   Usage: TF = SOFAcompare(A, B)
%
%   SOFAcompare compares A and B and
%   returns logical 1 (true) if they are identical and 0 (false)
%   if a difference has been found. 
%
%   [TF,reason,where] = SOFAcompare(..) provides in reason the reason 
%   for a difference and describes in where in which metadata the difference arose. 
%
%   SOFAcompare(A, B, 'ignoreDate') ignores the global attributes
%   DateCreated and DateModified. 
%
%   Currenty the functionality is limited to the comparison og 
%   attributes only.

% #Author: Piotr Majdak
% #Author: Michael Mihocic: header documentation updated (20.10.2021)
%
% SOFA Toolbox - function SOFAcompare
% Copyright (C) Acoustics Research Institute - Austrian Academy of Sciences
% Licensed under the EUPL, Version 1.2 or - as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "License")
% You may not use this work except in compliance with the License.
% You may obtain a copy of the License at: https://joinup.ec.europa.eu/software/page/eupl
% Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing  permissions and limitations under the License. 

definput.flags.type={'all','ignoreDate'};
[flags,~]=SOFAarghelper({},definput,varargin);


tf=1;
reason='';
where='';

%   % check if the size is equal
% o1=whos('Obj1');
% o2=whos('Obj2');
% if o1.bytes ~= o2.bytes, tf=0; reason='Different size'; return; end

  % get the field names
Xf=fieldnames(Obj1);

  % ignore DateCreated and DateModified?
if flags.do_ignoreDate
  Xf=fieldnames(rmfield(Obj1,{'GLOBAL_DateCreated','GLOBAL_DateModified'}));
end
  
  % check if we have the same fields in Obj2 as in Obj1
for ii=1:length(Xf)
  if ~isfield(Obj2,Xf{ii}), tf=0; reason='Field missing in B'; where=Xf{ii}; return; end
end

  % check if we have the same content of attributes in Obj2 as in Obj1
for ii=1:length(Xf)
  if isempty(strfind(Xf{ii},'_')), continue; end
  if ~strcmp(Obj1.(Xf{ii}),Obj2.(Xf{ii})), tf=0; reason='Field not equal'; where=Xf{ii}; return; end
end

