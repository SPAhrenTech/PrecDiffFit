classdef Sim<handle
    properties
        name;
        region;
        diff;
        H_refL=             [];
        intensL=            [];
        doPlot=             false;
        gMax=               1e4;
        sgMax=              0.1;
        invwMin=            0.02;
        useAllFitBeams=     true;
        lu_plotL=           [];
        eigenL=             {};
        dyn_luM=            [];
        dynH_refL=          [];
        usePert=            false;
        needsDiag=          true;
        useReflection=      false;%Calculate half precession and reflect
        beamPairL=          [];
        
        useCoher=           false;
        useIncoher=         false;
        dynamical=          false;
        fig;
        img;
        figNum
    end
    
    methods
        function sim=Sim(name)
            sim.name=name;
            %this.blochInfoL=cell_array;
        end
        
        function sim=createEigen(sim,nSteps,nBeams)
            for ii=1:nSteps
                sim.eigenL{ii,1}=Eigen(nBeams);
            end
        end
        
        %
        function sim=expandBeamL(sim,const)
            diff=sim.diff;
            region=sim.region;
            nZone=size(diff.N_ZoneL,2);
            %Expand sim beam list
            nDomain=size(region.domainL,1);
            for iDomain=1:nDomain
                domain=region.domainL{iDomain};
                hPrimMax=ceil(sim.gMax*norm(domain.A_primM'*const.nX));
                kPrimMax=ceil(sim.gMax*norm(domain.A_primM'*const.nY));
                lPrimMax=ceil(sim.gMax*norm(domain.A_primM'*const.nZ));
                
                % testM=region.A_refM*domain.B_primM'
                nBeamSim=size(sim.H_refL,1);
                %Add beams not in data list
                for hPrim=-hPrimMax:hPrimMax
                    for kPrim=-kPrimMax:kPrimMax
                        for lPrim=-lPrimMax:lPrimMax
                            H_primV=[hPrim kPrim lPrim]';
                            H_refV=region.A_refM*domain.B_primM'*H_primV;
                            gV=region.B_refM'*H_refV;
                            
                            useBeam=1;
                            for iBeamSim=1:nBeamSim%check if already included
                                Hp_refV=sim.H_refL(iBeamSim,:)';
                                gpV=region.B_refM'*Hp_refV;
                                dgV=gpV-gV;
                                if norm(dgV)<1e-7
                                    useBeam=0;
                                    break;
                                end
                            end
                            
                            if ~useBeam
                                continue
                            end
                            
                            if norm(gV)>sim.gMax+1e-8%truncation issue
                                continue;
                            end
                            
                            %{
                            sg=getExcErr(gV,diff.kV,const.nZ);
                            if(abs(sg)>sim.sgMax)
                                iXsig=diff.lambda*getUg(const,sim,H_refV);
                                invw=norm(iXsig/sg);
                                if invw<sim.invwMin
                                    continue
                                end
                            end
                            %}
                            
                            %Check if beam is in included zones
                            useBeam=0;
                            for iZone=1:nZone
                                if(abs(dot(H_refV,region.uvwZoneV)-diff.N_ZoneL(iZone))<1e-7)
                                    useBeam=1;
                                end
                            end
                            
                            %Add the beam to the lists
                            if useBeam
                                nBeamSim=nBeamSim+1;
                                sim.H_refL(nBeamSim,:)=H_refV;
                            end
                            
                        end
                    end
                end
            end
            
            %Make sure sim beam list contains [0 0 0]; (it should)
            simContains000=false;
            for iBeamSim=1:nBeamSim
                H_refV=sim.H_refL(iBeamSim,:)';
                gV=region.B_refM'*H_refV;
                if norm(gV)<1e-7
                    simContains000=true;
                    break;
                end
            end
            
            if ~simContains000
                sim.H_refL=[[0 0 0];sim.H_refL];
                nBeamSim=nBeamSim+1;
            end
            
            %%Sort sim beams by length of g
            gL=region.B_refM'*sim.H_refL';
            len_gL=vecnorm(gL)';
            [len_gL,sortL]=sort(len_gL);
            
            Hp_refL=sim.H_refL(sortL(:,1),:);
            sim.H_refL=Hp_refL;
        end
        
         %Find Friedel pairs (assumes ZOLZ)
        function sim=findPairs(sim)        
            nBeam=size(sim.H_refL,1);
            sim.beamPairL=[1:nBeam]';
            for ii=1:nBeam
                H_refV1=sim.H_refL(ii,:);
                g1V=H_refV1*sim.region.B_refM;

                for jj=ii+1:nBeam
                    H_refV2=sim.H_refL(jj,:);
                    g2V=H_refV2*sim.region.B_refM;
                    
                    gAvgV=g1V+g2V;
                    if norm(gAvgV)<1e-7 %%Associate pairs                        
                        sim.beamPairL(ii,1)=jj;
                        sim.beamPairL(jj,1)=ii;
                     end
                end
            end
        end   

        
        function sim=getBeamL(sim,const,fit)
            
            expt=fit.expt;
            region=sim.region;
            diff=sim.diff;
            
            nZone=size(diff.N_ZoneL,1);
            %Create sim beam list from compatible expt beams
            sim.H_refL=[];
            nDomain=size(region.domainL,1);
            nBeamSim=0;
            nBeamExpt=size(expt.intensL,1);
            for iBeamExpt=1:nBeamExpt
                H_refV=expt.H_refL(iBeamExpt,:)';
                gV=region.B_refM'*H_refV;
                
               addBeam=false;              
               for iDomain=1:nDomain
                    domain=region.domainL{iDomain};
                    H_primV=domain.A_primM*gV;
                    dH_primV=H_primV-round(H_primV);
                    if norm(dH_primV)<1e-7
                        addBeam=true;                        
                    end
               end
               if ~addBeam
                    continue;
                end               
                addBeam=false;

                
                len_g=norm(gV);
                if sim.useAllFitBeams
                    if (fit.gMin<=len_g)&&(len_g<=fit.gMax)
                        addBeam=true;
                     end
                end
                                 
                 if ~addBeam
                    if len_g<=sim.gMax
                        for iZone=1:nZone
                            if(abs(dot(H_refV,region.uvwZoneV)-diff.N_ZoneL(iZone))<1e-7)
                            	addBeam=true;
                            end
                        end
                     end
                 end
                 
                if addBeam
                   nBeamSim=nBeamSim+1;
                   sim.H_refL(nBeamSim,:)=H_refV';
                end

            end
            
            nPhase=size(region.domainL(),1);
            for iPhase=1:nPhase
                phase=region.domainL{iPhase}.phase;
                phase.adjust();
            end
            
            sim.expandBeamL(const);
            sim.findPairs();
            %sim.createBlochInfo(diff.nSteps,nBeamSim);
            sim.createEigen(diff.nSteps,nBeamSim);
        end
        
        %Structure factors, separated into Hermitian and non-Hermitian
        %terms as [UgH UgNH]
        function Ug=getUg(sim,const,H_ref)
            
            region=sim.region;
            basis=region.basis;
            diff=sim.diff;
            gV=basis.B_refM'*H_ref;
            s=norm(gV)/2;
            Ug_V=[0 0];%Hermitian and non-hermitian combined
            totVol=0;
            nDomain=size(region.domainL,1);
            for iDomain=1:nDomain
                domain=region.domainL{iDomain};
                phase=domain.phase;
                H_prim=phase.A_primM*gV;
                dH_prim=H_prim-round(H_prim);
                
                totVol=totVol+domain.relVol;
                
                if norm(dH_prim)>1e-7
                    continue
                end
                
                nElem=size(phase.elemL,2);
                Fg=0;
                for iElem=1:nElem
                    fe=phase.occL(iElem)*get_fe(const,diff,phase.elemL(iElem),s)...
                        *[1 1i*phase.absorp];
                    Xm=phase.adjCoordL(:,phase.siteL(iElem));
                    dV=basis.A_refM'*Xm;
                    dFg=fe*exp(-2*pi*1i*(dot(gV,dV)));
                    Fg=Fg+dFg;
                end
                vPrim=det(phase.A_primM);
                dUg=1/pi/vPrim*exp(-phase.B_T*norm(gV)^2/4)*Fg;
                Ug_V=Ug_V+domain.relVol*dUg;
            end
            if totVol==0
                Ug=0;
            else
                Ug=Ug_V/totVol;
            end
        end
        
        %Structure factors: just the Hermitian part
        function Ug=getUgH(sim,const,H_ref)
            UgF=sim.getUg(const,H_ref);
            Ug=UgF(1,:);
        end
        
        %Fill out list of referenced structure factor list
        function sim=getDynH_refL(sim)
            
            nSimRef=size(sim.H_refL,1);
            %nBeam=size(beamList,1);
            sim.dynH_refL=[];
            sim.dyn_luM=zeros(nSimRef,nSimRef);
            nRef=0;
            for ii=1:nSimRef
                H_refVi=sim.H_refL(ii,:)';
                for ji=1:nSimRef
                    H_refVj=sim.H_refL(ji,:)';
                    dH_refV=H_refVi-H_refVj;
                    
                    gV=sim.region.B_refM'*dH_refV;
                    addRef=true;
                    for iRef=1:nRef%check if already included
                        Hp_refV=sim.dynH_refL(iRef,:)';
                        gpV=sim.region.B_refM'*Hp_refV;
                        dgV=gpV-gV;
                        if norm(dgV)<1e-7
                            addRef=false;
                            ijRef=iRef;
                            break;
                        end
                    end
                    
                    if(addRef)
                        nRef=nRef+1;
                        sim.dynH_refL(nRef,:)=dH_refV;
                        ijRef=nRef;
                    end
                    sim.dyn_luM(ii,ji)=ijRef;
                end
            end
        end
        
      %Fill out list of referenced structure factor list
        function sim=getDynH_refL2(sim)
            
            nSimRef=size(sim.H_refL,1);
            sim.dynH_refL=[];
            full_dgRefL=[];
            sim.dyn_luM=zeros(nSimRef,nSimRef);
            nRef=0;
            for ii=1:nSimRef
                H_refVi=sim.H_refL(ii,:);
                for ji=1:nSimRef
                    H_refVj=sim.H_refL(ji,:);
                    dH_refV=H_refVi-H_refVj;
                    dgRefV=dH_refV*sim.region.B_refM;
                    exists=false;
                    if nRef>0
                       dgRefL=repmat(dgRefV,nRef,1);
                       [~,list]=ismembertol(dgRefL,full_dgRefL,1e-7,'ByRows',true);
                       ijRefL=find(list,1,'first');%check if already included
                       if size(ijRefL,1)>0
                           ijRef=ijRefL(1);  
                           exists=true;
                        end
                    end
                    if ~exists
                        nRef=nRef+1;
                        sim.dynH_refL(nRef,:)=dH_refV;
                        full_dgRefL(nRef,:)=dgRefV;
                        ijRef=nRef;
                    end
                    sim.dyn_luM(ii,ji)=ijRef;
                end
            end
       end
        %Fill out Hermitian matrix iXsiMat with inverse extinction distances
        function iXsiM=get_iXsiM(sim,const)
            nDyn=size(sim.dynH_refL,1);
            iXsiL=zeros(nDyn,2);
            for ii=1:nDyn
                dH_ref=sim.dynH_refL(ii,:)';
                Ug=sim.getUg(const,dH_ref);
                iXsiL(ii)=sim.diff.lambda*(Ug(1,1)+Ug(1,2));
            end
            nBeam=size(sim.H_refL,1);
            iXsi0M=zeros(nBeam,nBeam);
            iXsi0M=iXsiL(sim.dyn_luM);
            %for ii=1:nBeam
            % for ij=1:nBeam
            %iXsi0M(ii,ij)=iXsi0L(sim.dyn_luM(ii,ij),1);
            %end
            %end
            %[ilu,~]=ngrid(1:nBeam,1:nBeam);
            %iXsi0M(sub2ind(size(iXsi0M),ilu,sim.dyn_luM))=2;
            %iXsi0L(sim.dyn_luM(:,:),1);
            % iXsi0M=iXsi0L(sim.dyn_luM(:,:),:);
            %Separate Hermitian and Non-hermitian matrices
            iXsiM=zeros(nBeam,nBeam,2);
            iXsiM(:,:,1)=(iXsi0M+iXsi0M')/2;
            iXsiM(:,:,2)=(iXsi0M-iXsi0M')/2/1i;
        end
        
        
         
        %Fill out Hermitian matrix iXsiMat with inverse extinction distances
        function UgL=getUgL(sim,const)
            nBeam=size(sim.H_refL,1);
            UgL=zeros(nBeam,2);
            for ii=1:nBeam
                H_ref=sim.H_refL(ii,:)';
                U=sim.getUg(const,H_ref);
                UgL(ii,:)=U;
            end
        end
        
        %Appropriate parameters should already be set
        function sim=prepare(sim,const)
            if ~sim.dynamical
                return
            end
            nBeams=size(sim.H_refL,1);
            diff=sim.diff;
            iXsiM=sim.getDynH_refL().get_iXsiM(const);
            first=true;
            for kk=1:diff.nSteps
               if sim.useReflection
                   angPrec=pi*kk/diff.nSteps+pi/2;
               else
                    angPrec=2*pi*kk/diff.nSteps+pi/2;
               end
               angPrec=2*pi*kk/diff.nSteps+pi/2;
                phiX=           diff.phi*cos(angPrec);
                phiY=           diff.phi*sin(angPrec);
                nkX=            sin(phiX/1000);
                nkY=            sin(phiY/1000);
                nkZ=            sqrt(1-nkX^2-nkY^2);
                nkV=            nkX*const.nX+nkY*const.nY+nkZ*const.nZ;%beam direction
                diff.kV=        nkV/diff.lambda;
                
                sim.eigenL{kk,1}=Eigen(nBeams);
                sim.eigenL{kk,1}.prepare(sim,iXsiM);
                %sim.eigenL{kk}=etemp;
            end
        end
        
        %
        function summarize(sim)
            nBeams=0;
            if sim.dynamical               
                nStrong=0;nPert=0;nNegl=0;            
                for kk=1:sim.diff.nSteps
                    nStrong=nStrong+size(sim.eigenL{kk}.strongL,1);;
                    nPert=nPert+size(sim.eigenL{kk}.pertL,1);
                    nNegl=nNegl+size(sim.eigenL{kk}.negL,1);
                end
               nBeams=nStrong+nPert+nNegl;           
               fStrong=nStrong/nBeams;
               fPert=+nPert/nBeams;
               fNegl=nNegl/nBeams;
                sSumm=compose("total beam calculations: %6.0f",[nBeams]);sFormat="\n%s";
                
                sSummP=compose("strong: %5.2f%%, perturbative: %5.2f%%, neglected: %5.2f%%",100*[fStrong fPert fNegl]);              
                sSumm=[sSumm sSummP];
                sFormat=sFormat+"\n%s";

                nU=size(sim.dynH_refL,1);
                sSummP=compose("structure factors: %6.0f",nU);              
                sSumm=[sSumm sSummP];
                sFormat=sFormat+"\n%s";

               sprintf(sFormat,sSumm)          
            else
               nBeams=size(sim.H_refL,1);
               sSumm=compose("total beam calculations: %6.0f",[nBeams]);sFormat="\n%s";           
               sprintf(sFormat,sSumm)          
           end
          end
        
        % simList contains indices of beams to include
        function psiM=getPsiM(sim,eigen)
            [nBeam,nBloch]=size(eigen.C_M);
            [col,row]=meshgrid(1:nBloch,1:nBeam);
            dExpM=eigen.gammaL(col)-eigen.sgL(row);
            phaseFacM=exp(2*pi*1i*dExpM*sim.region.thick);
            epsM=repmat(conj(eigen.C_M(1,:)),nBeam,1);
            psiM=eigen.C_M.*phaseFacM.*epsM;%.*coherM;
        end
        
        % simList contains indices of beams to include
        function psiM=getPsiCoherM(sim,eigen)
            [nBeam,nBloch]=size(eigen.C_M);
            [col,row]=meshgrid(1:nBloch,1:nBeam);
            
            dExpM=eigen.gammaL(col)-eigen.sgL(row);
            phaseFacM=exp(2*pi*1i*dExpM*sim.region.thick);
            coherM=exp(-pi^2*(dExpM*sim.region.dCoherThick).^2);
            phaseFacM=coherM.*phaseFacM;
            epsM=repmat(conj(eigen.C_M(1,:)),nBeam,1);
            psiM=eigen.C_M.*phaseFacM.*epsM;%.*coherM;
        end
        
        % simList contains indices of beams to include
        function I_L=getIntensL(sim,psiM)
            I_L=abs(diag(psiM*psiM'));
        end
        
        % simList contains indices of beams to include
        function I_L=getIntensIncoherL(sim,psiM,eigen)
            nBloch=size(psiM,2);
            col=meshgrid(1:nBloch,1:nBloch);
            gammaM=eigen.gammaL(col);
            dGammaM=gammaM-gammaM.';
            incoherM=exp(-pi^2*(dGammaM.*sim.region.dIncoherThick).^2);
            I_L=abs(diag(psiM*incoherM*psiM'));
        end
        
        % simList contains indices of beams to include
        function psigL=getPsigL(sim,eigen)
            gammaM=diag(exp(2*pi*1i*eigen.gammaL*sim.region.thick));
            epsV=eigen.psiM(1,:)';
            psigL=getPsiCoherM(sim,eigen)*gammaM*epsV;
        end
        
        % simList contains indices of beams to include
        function I_L=getIntensKinL(sim,const)
            UgL=sim.getUgL(const);
            I_L=abs(UgL(:,1)).^2;
        end
        
               
        %
        function sim=getIntensPrecL(sim,const)
            diff=sim.diff;
            
            if ~sim.dynamical
                if sim.needsDiag
                    sim.intensL=sim.getIntensKinL(const);
                end
                return
            end
            if sim.needsDiag                
                iXsiM=sim.get_iXsiM(const);
                for kk=1:diff.nSteps
                    
                    if sim.useReflection
                        angPrec=pi*kk/diff.nSteps+pi/2;
                    else
                        angPrec=2*pi*kk/diff.nSteps+pi/2;
                    end
                    phiX=           diff.phi*cos(angPrec);
                    phiY=           diff.phi*sin(angPrec);
                    nkX=            sin(phiX/1000);
                    nkY=            sin(phiY/1000);
                    nkZ=            sqrt(1-nkX^2-nkY^2);
                    nkV=            nkX*const.nX+nkY*const.nY+nkZ*const.nZ;%beam direction
                    diff.kV=        nkV/diff.lambda;
                    
                    sim.eigenL{kk,1}=sim.eigenL{kk,1}.calculate(const,sim,iXsiM);
                end
            end
            sim.needsDiag=false;
            
            first=true;
            for kk=1:diff.nSteps
                
                if sim.useCoher
                    psiM=sim.getPsiCoherM(sim.eigenL{kk,1});
                else
                    psiM=sim.getPsiM(sim.eigenL{kk,1});
                end
                
                if sim.useIncoher
                    intensL=sim.getIntensIncoherL(psiM,sim.eigenL{kk,1});
                else
                    intensL=sim.getIntensL(psiM);
                end
                
                %psigL=sim.getPsigL(const,sim.eigenL{kk,1});
                %intensL=conj(psigL).*psigL;
                
                if first
                    %sim.eigenL{kk}.psiM
                    sim.intensL=intensL;
                else
                    switch diff.mode
                        case 'max'%maximum
                            sim.intensL(intensL>sim.intensL)=intensL(intensL>sim.intensL);
                            
                        case 'int'%integrated
                            sim.intensL=sim.intensL+intensL;
                    end                    
                end
                
                if sim.useReflection
                    intensLp=intensL(sim.beamPairL(:,1),1);
                    switch diff.mode
                        case 'max'%maximum
                            sim.intensL(intensLp>sim.intensL)=intensLp(intensLp>sim.intensL);
                            
                        case 'int'%integrated
                            sim.intensL=sim.intensL+intensLp;
                    end
                end
                first=false;
            end
            
        end
        
        %Create plot list
        function sim=createPlotList(sim,luL,plotFitBeamsOnly,gPlotMax)
            nBeam=size(sim.H_refL,1);
            sim.lu_plotL=[];
            nPlotBeam=0;
            for iBeam=1:nBeam
                H2_refV=sim.H_refL(iBeam,:)';
                g2V=sim.region.B_refM'*H2_refV;
                doPlotBeam=true;
                if norm(g2V)>gPlotMax
                    doPlotBeam=false;
                end
                if (plotFitBeamsOnly)&&(~ismember(iBeam,luL))
                    doPlotBeam=false;
                    %sim.lu_plotL=luL;
                end
                if doPlotBeam
                    nPlotBeam=nPlotBeam+1;
                    sim.lu_plotL(nPlotBeam,1)=iBeam;
                end
            end
        end
        
        function sim=plot(sim,const)
            if isempty(sim.fig)
                sim.fig=figure('Name',sim.name,'NumberTitle','off');
                sim.figNum=get(gcf,'Number');
            end
            if sim.diff.doRealPlot
                sim.img=sim.diff.realPlot(const,sim.region,sim.intensL,sim.H_refL,sim.lu_plotL);
            else
                sim.img=sim.diff.plot(const,sim.region,sim.intensL,sim.H_refL,sim.lu_plotL);
            end 
            showFig=false;
            if isvalid(gcf)
                frontFigNum=get(gcf,'Number');
                if frontFigNum~=sim.figNum
                    showFig=true;
                end
            end
            if showFig
                figure(sim.figNum);
            end
            imshow(sim.img,'InitialMagnification',100);
        end
        
    end
end