function [thickResL,I_resL]=thicknessScan(const,basis,sim,peakInfo,scanInfo)

    diff=sim.diff;
    region=sim.region;
 
    nPhase=size(region.domainL(),1);
    for iPhase=1:nPhase        
        phase=region.domainL{iPhase}.phase;
        phase.adjust();   
    end
    
    UgBeam=getUg(const,sim,peakInfo.H_refBragg)
    XsiBeam=1/(sim.diff.lambda*UgBeam)
    
   %Get beams
    sim.H_refL=[];
    sim.H_refL=[0 0 0;peakInfo.H_refBragg'];
    expandSimBeamL(const,sim,region);
    nBeamSim=size(sim.H_refL,1);
    lu_simL=[1:nBeamSim]';

    %Find Bragg beam
    nBeamSim=size(sim.H_refL,1);
    for iBeamSim=1:nBeamSim%check if already included
        Hp_refV=sim.H_refL(iBeamSim,:)';
        dHV=Hp_refV-peakInfo.H_refBragg;
        if norm(dHV)<1e-7
            peakInfo.iBragg=iBeamSim;
            break;
        end
    end
    
    UgBeam=getUg(const,sim,peakInfo.H_refBragg)
    XsiBeam=1/(sim.diff.lambda*UgBeam)

    %Dynamical part 
    sim.createBlochInfo(1,nBeamSim);
    iXsiM=getiXsiM(const,sim);
    sim.blochInfoL{1}=sim.blochInfoL{1}.getBlochInfo(const,sim,iXsiM);

    IL=[1:scanInfo.nThick];
    thickL=[1:scanInfo.nThick];

    thickScan=exp(log(scanInfo.thickMax/scanInfo.thickMin)/(scanInfo.nThick-1));
    for iThick=1:scanInfo.nThick

        thick=scanInfo.thickMin*thickScan^(iThick-1); %nm%

        %thickness
        region.thick=              thick;
        psigL=getPsigL(const,sim,sim.blochInfoL{1});
        sim.intensL=conj(psigL).*psigL;
        I0=sim.intensL(peakInfo.iBragg);

        thickL(iThick)=thick;
        IL(iThick)=I0;
    end

    hold on;
    p=plot(thickL,IL);
    set(p,'linewidth',[4.0])
    set(gca,'Fontsize',[12])
    xlabel('thickness (nm)')
    ylabel('Intensity')

    grid on

    thickResL=thickL';
    I_resL=IL';

end
