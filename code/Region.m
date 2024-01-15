classdef Region<handle
    properties
        name='';
        domainL=[];
        A_refM=[];
        B_refM=[];
        orientM=                [1 0 0;0 1 0;0 0 1];
        basis;
        uvwZoneV=               [0 0 1]';
        hklXV=                  [1 -1 0]';
        relVolTot;
        etaV=[];
        volFracV=[];
        thick=                  0;%thickness
        
        coherThickRat=              0.1;%Coherent fraction of thickness variation
        dCoherThick=            0%coherent thickness

        incoherThickRat=            1;%Fractional variation in thickness
    	dIncoherThick=         0;%thickness variation
        
    end
    
    properties (Dependent)
     % phaseL;
    end
    
    methods
        function h=Region(name,basis)
            if nargin > 0
                h.name=name;
                h.basis=basis;
                h.A_refM=basis.A_refM;
            end
        end
                
        %
        function L=phaseL(region)
            nDomain=size(region.domainL,1);
            L={};
            for iDomain=1:nDomain
                newPhase=true;
                domain=region.domainL{iDomain};
                phase=domain.phase;
                nPhase=size(L,1);
                for iPhase=1:nPhase
                    if phase==L{iPhase}
                        newPhase=false;
                        break;
                    end
                end
                if newPhase
                    L{nPhase+1,1}=phase;
                end
            end
        end
        
      function region=addDomain(region,domain,relVol)
           nDomain=size(region.domainL,1);
           region.domainL{nDomain+1,1}=domain;
           domain.region=region;
           domain.phase.basis=region.basis;
           domain.relVol=relVol;
      end
      
      function region=adjust(region)
         region.dCoherThick=region.coherThickRat*region.thick;
        region.dIncoherThick=region.incoherThickRat*region.thick;
        end
       
        function region=alignVectors(region,oldV,newV)
            %Orient reference
            RM=getRotationAlignVectors(oldV,newV);
            region.A_refM=region.A_refM*RM';
            region.B_refM=inv(region.A_refM)';

            %Ordnet domains
            nDomain=size(region.domainL,1);
            for iDomain=1:nDomain
                domain=region.domainL{iDomain};
                domain.A_primM=domain.A_primM*RM';
                domain.B_primM=inv(domain.A_primM)';
            end

            %Save orientation matrix
            region.orientM=region.orientM*RM';

        end

       function region=orient(region,const,diff)
            %Rotate zone axis to foil normal (z)
            rZoneV=region.A_refM'*region.uvwZoneV;
            region.alignVectors(rZoneV,diff.nFoilV);

            %Rotate hklX to x axis
            gX_V=region.B_refM'*region.hklXV;
            region.alignVectors(gX_V,const.nX);

            region.B_refM=inv(region.A_refM)';

            %Orient domains
            nDomain=size(region.domainL,1);
            for iDomain=1:nDomain
                domain=region.domainL{iDomain};
                domain.B_primM=inv(domain.A_primM)';
            end

       end

       %
       function region=getVolFracV(region)
            nDomainP=size(region.domainL,1);
            region.etaV=zeros(nDomainP,1);
            region.volFracV=zeros(nDomainP,1);
            for iDomain=1:nDomainP
                domain=region.domainL{iDomain};
                region.volFracV(iDomain,1)=domain.relVol;
                phase=domain.phase;
                if strcmp(phase.type,'ordered')
                    region.etaV(iDomain,1)=phase.eta;
                end
            end
            region.relVolTot=sum(region.volFracV);
            region.volFracV=region.volFracV/region.relVolTot;       
       end
        
       function eta=getEta(region)
            eta=0;
            nDomainP=size(region.domainL,1);
            for iDomain=1:nDomainP
                domain=region.domainL;
                eta=eta+region.etaV(iDomain)*region.volFracV(iDomain);
            end
       end
    end
end
