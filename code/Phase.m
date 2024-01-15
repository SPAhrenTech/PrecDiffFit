classdef Phase<handle
   properties
	name='';
	A_primM=[];%Not rotated
    basis;
	elemL={};
	siteL=[];
	occL=[];
    coordL=[];
    adjCoordL=[];
	type='';
	B_T=         	0.0020;%(1/nm^2) Debye-Waller factor. Currently just a guess
    absorp=       	0.1;%ratio of non-hermitian/hermitian structure factors

   end
	methods
        function h=Phase(nameP)
            if nargin>0
                h.name = nameP;
            end
        end
        
       function adjust(phase)
            phase.adjCoordL=phase.coordL;
       end
	end
end
