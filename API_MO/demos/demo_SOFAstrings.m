X=SOFAload([SOFAdbPath '\SOFA\ARI_NH2_hrtf_M_dtf 256.sofa']);
str={};
for ii=1:X.API.M
  str{ii,1}=['NH' num2str(round(rand(1,1)*10000))];
end
X2=SOFAaddVariable(X,'Test','MS',str);    
X2=SOFAupdateDimensions(X2);
SOFAsave('test.sofa',X2);

%%
N=SOFAload('test.sofa');