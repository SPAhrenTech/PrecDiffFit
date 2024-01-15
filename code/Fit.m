classdef Fit<handle
    properties
        expt;
        sim;
        region;
        diff;
        luL=[];%col 1=expt;col 2=sim
        plotFitBeamsOnly=      false;%whether to only plot fitted beams
        gMax=                   1e3;%(1/nm)
        gMin=                   0;%(1/nm)
        fitF=                   false;%fit structure factors: true, intensities: false
        H_specRefL=             [];%Specific list of beams to include
        fitPairsOnly=           false;

        valL=                   [];%expt;sim%sig;w      
    end
    
    methods
        function fit=Fit(exptP,simP,regionP,diffP)
            if nargin>0
                fit.expt=exptP;
                fit.sim=simP;
                fit.region=regionP;
                fit.diff=diffP;
                fit.sim.diff=diffP;
                fit.sim.region=regionP;
                fit.expt.diff=diffP;
                fit.expt.region=regionP;
                fit.luL=[];
            end
        end
        
        %Read a text file containg a text hearder line followed by
        %four columns: h k l Intens
        function fit=readExptBeamL(fit)
            fit.expt.readBeamL(fit.region);
        end
        
        function fit=getSimBeamL(fit,const)
            fit.sim.getBeamL(const,fit);
        end
        
        %
        function fit=scaleIntens(fit,const)           
           if fit.fitF%structure factors
              fit.valL(:,2)=sqrt(fit.sim.intensL(fit.luL(:,2),1)); 
           else
              fit.valL(:,2)=fit.sim.intensL(fit.luL(:,2),1);            
           end
           
           %Find scaling ratio
           numer=sum(fit.valL(:,4).*fit.valL(:,1).*fit.valL(:,2)./fit.valL(:,3).^2);
           denom=sum(fit.valL(:,4).*fit.valL(:,2).^2./fit.valL(:,3).^2);
           ratio=numer/denom;
           fit.valL(:,2)=ratio*fit.valL(:,2);
           
           if fit.fitF
             fit.sim.intensL=ratio^2*fit.sim.intensL;            
           else
              fit.sim.intensL=ratio*fit.sim.intensL;            
           end 
        end
        
        %Create expt<->sim look-up lists
        function fit=createSpecBeamLists(fit)
            expt=fit.expt;
            sim=fit.sim;
            nBeamSpec=size(fit.H_specRefL,1);
            nBeamExpt=size(fit.expt.H_refL,1);
            nBeamSim=size(fit.sim.H_refL,1);
            nBeamFit=size(fit.luL,1);
            for iBeamSpec=1:nBeamSpec
                H_specRefV=fit.H_specRefL(iBeamSpec,:)';
                g_specV=fit.region.B_refM'*H_specRefV;
                
                %Check if in expt
                addBeam=false;
                for iBeamExpt=1:nBeamExpt
                    H_exptRefV=expt.H_refL(iBeamExpt,:)';
                    g_exptV=fit.region.B_refM'*H_exptRefV;
                    dgV=g_specV-g_exptV;
                    if norm(dgV)<1e-7
                        addBeam=true;
                        break;
                    end
                end
                
                %Check if in sim
                if addBeam
                    addBeam=false;
                    for iBeamSim=1:nBeamSim
                        H_simRefV=sim.H_refL(iBeamSim,:)';
                        g_simV=fit.region.B_refM'*H_simRefV;
                        dgV=g_specV-g_simV;
                        if norm(dgV)<1e-7
                            addBeam=true;
                            break;
                        end
                    end
                end
                
                if addBeam
                    nBeamFit=size(fit.luL,1);
                    for iBeamFit=1:nBeamFit
                        if (iBeamExpt==fit.luL(iBeamFit,1))&&(iBeamSim==fit.luL(iBeamSim,1))
                            addBeam=false;
                            break;
                        end
                    end
                end
                
                if addBeam
                    nBeamFit=size(fit.luL,1);
                    nBeamFit=nBeamFit+1;
                    fit.luL(nBeamFit,:)=[iBeamExpt iBeamSim];
                end
                
            end
        end
        
        %Create expt<->sim look-up lists
        function fit=createBeamLists(fit)
            expt=fit.expt;
            sim=fit.sim;
            nBeamExpt=size(fit.expt.H_refL,1);
            nBeamSim=size(fit.sim.H_refL,1);
            nBeamFit=size(fit.luL,1);
            for iBeamExpt=1:nBeamExpt
                addBeam=false;
                if fit.fitPairsOnly&&(~expt.beamPairL(iBeamExpt,1))
                    continue
                end
                H_exptRefV=expt.H_refL(iBeamExpt,:)';
                g_exptV=fit.region.B_refM'*H_exptRefV;
                len_g_exptV=norm(g_exptV);
                if (len_g_exptV<=fit.gMax)&&(len_g_exptV>=fit.gMin)
                    for iBeamSim=1:nBeamSim
                        H_simRefV=sim.H_refL(iBeamSim,:)';
                        g_simV=fit.region.B_refM'*H_simRefV;
                        dgV=g_simV-g_exptV;
                        if norm(dgV)<1e-7
                          addBeam=true;                               
                         break;
                        end
                    end
                end
                
                if addBeam
                    nBeamFit=size(fit.luL,1);
                    for iBeamFit=1:nBeamFit
                        if (iBeamExpt==fit.luL(iBeamFit,1))&&(iBeamSim==fit.luL(iBeamSim,1))
                            addBeam=false;
                            break;
                        end
                    end
                end
                
                if addBeam
                    nBeamFit=nBeamFit+1;
                    fit.luL(nBeamFit,:)=[iBeamExpt iBeamSim];
                end
            end
            
            fit.createSpecBeamLists();
            nBeamFit=size(fit.luL,1);
            
            fit.valL=zeros(nBeamFit,7);%All the fit values are stored here.            
            if nBeamFit>0
                if(fit.fitF)
                   fit.valL(:,1)=sqrt(expt.intensL(fit.luL(:,1)));
                   fit.valL(:,3)=expt.sigIL(fit.luL(:,1))/2./fit.valL(:,1);
                else
                    fit.valL(:,1)=expt.intensL(fit.luL(:,1));
                    fit.valL(:,3)=expt.sigIL(fit.luL(:,1));
                end
                fit.valL(:,4)=ones(nBeamFit,1);%weight
                fit.valL(:,5:7)=sim.H_refL(fit.luL(:,2),:);
            else
                compose("No beams to fit.")
            end
        end
               
        %Adjust weights of beams that would arise from compPhase
        %by factor wFac
        function fit=adjustWeightL_byPhase(fit,compPhase,wFac,exclusive)
            expt=fit.expt;
            
            nBeamFit=size(fit.luL,1);
            for iBeamFit=1:nBeamFit
                H_fitRefV=expt.H_refL(fit.luL(iBeamFit,1),:)';
                g_fitRefV=fit.region.basis.B_refM'*H_fitRefV;
                H_primV=compPhase.A_primM*g_fitRefV;
                dH_primV=H_primV-round(H_primV);
                doChange=(exclusive==(norm(dH_primV)<1e-7));
                if doChange
                    fit.valL(iBeamFit,4)=wFac;
                end
            end
        end
        
          %Adjust weights of beams that would arise from compPhase
        %by factor wFac
        function fit=adjustWeightL_byRefL(fit,H_compRefL,wFac,exclusive)
            expt=fit.expt;
            nBeamComp=size(H_compRefL,1);
            nBeamFit=size(fit.luL,1);
            for iBeamFit=1:nBeamFit
                H_fitRefV=expt.H_refL(fit.luL(iBeamFit,1),:)';
                g_fitRefV=fit.region.basis.B_refM'*H_fitRefV;
                for iBeamComp=1:nBeamComp
                    H_compRefV=H_compRefL(iBeamComp,:)';
                    g_compRefV=fit.region.basis.B_refM'*H_compRefV;
                    dgV=g_compRefV-g_fitRefV;
                    doChange=(exclusive==(norm(dgV)<1e-7));
                    if doChange
                        fit.valL(iBeamFit,4)=wFac;
                    end
                end
            end
        end
      function X_V=getX_V(fit,const)
            fit.sim.getIntensPrecL(const);
            fit.scaleIntens(const);
            X_V=fit.valL(:,4).*(fit.valL(:,1)-fit.valL(:,2))./fit.valL(:,3);
            wAvg=sum(fit.valL(:,4))/size(fit.valL(:,4),1);
            X_V=X_V/wAvg;
        end
        
        function intensL=getSimIntensL(fit,const)
            fit.sim.getIntensPrecL(const);
            intensL=fit.sim.intensL;
        end
        
       
        %
        function fit=plotExpt(fit,const)
            expt=fit.expt;
            if expt.doPlot
                expt.plot(const);
            end
        end
        
        function fit=plotSim(fit,const)
            sim=fit.sim;
            if sim.doPlot
                if isempty(sim.fig)
                    sim.fig=figure('Name',sim.name,'NumberTitle','off');
                    sim.figNum=get(gcf,'Number');
                    expt=fit.expt;
                    if ~isempty(expt.fig)
                        simPos=expt.fig.Position;
                        simPos(1)=simPos(1)+simPos(3);
                        sim.fig.Position=simPos;
                    end
                end
                sim.plot(const);
            end
        end
        
    end
end
