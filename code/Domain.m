classdef Domain<handle
   properties
      region;
      phase;
      relVol;
      A_primM=[];
      B_primM=[];
   end
   methods
        function domain=Domain(phaseP)
            if nargin>0
                 domain.phase=phaseP;
                domain.A_primM=phaseP.A_primM;
              end
        end
   end
end
