function [tiltResL,I_resL]= rockingCurve(const,basis,sim,peakInfo,scanInfo)

    diff=sim.diff;
    region=sim.region;
    
    nPhase=size(region.domainL(),1);
    for iPhase=1:nPhase        
        phase=region.domainL{iPhase}.phase;
        phase.adjust();   
    end
    
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
    gBraggV=region.B_refM'*peakInfo.H_refBragg;
    peakInfo.tiltAxisV=cross(const.nZ,gBraggV);
    peakInfo.tiltAxisV=peakInfo.tiltAxisV/norm(peakInfo.tiltAxisV);

    %Dynamical part
    sim.createBlochInfo(1,nBeamSim);
    iXsiM=getiXsiM(const,sim);
 
    %rocking curve
    IL=[1:scanInfo.nTilt];
    tiltL=[1:scanInfo.nTilt];
    kV0=diff.kV;
    for iTilt=1:scanInfo.nTilt
        %Tilt beam slightly away from Bragg condition
        dTheta=scanInfo.tiltMin+(scanInfo.tiltMax-scanInfo.tiltMin)*(iTilt-1)/(scanInfo.nTilt-1);
        RM=getRotationAboutAxis(peakInfo.tiltAxisV,dTheta);%
        diff.kV=RM*kV0;
        sim.blochInfoL{1}=sim.blochInfoL{1}.getBlochInfo(const,sim,iXsiM);
        psigL=getPsigL(const,sim,sim.blochInfoL{1});
        sim.intensL=conj(psigL).*psigL;
        I0=sim.intensL(peakInfo.iBragg);

        tiltL(iTilt)=dTheta;
        IL(iTilt)=I0;
    end

    hold on;
    p=plot(tiltL,IL);
    set(p,'linewidth',[4.0])
    set(gca,'Fontsize',[12])
    xlabel('tilt (rad)')
    ylabel('Intensity')

    tiltResL=tiltL';
    I_resL=IL';

end

