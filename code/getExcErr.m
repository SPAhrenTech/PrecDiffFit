% --------------------------------------------------------------------------------
% Excitation error - P. Ahrenkiel (2018)
% g= RLV; k=incident beam wave vector; n= direction of s
% --------------------------------------------------------------------------------
%
function sg=getExcErr(gV,kV,nV)
   b=2*dot(kV+gV,nV);
   c=-(2*dot(kV,gV)+dot(gV,gV));
    
   sg=c/b;    
end   
 