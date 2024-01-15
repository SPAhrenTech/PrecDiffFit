classdef Proj<handle
    
    properties
        fitL={};
        const;
        basis;
        nIters=100;
        
        parL={};
        constrL={};
        varL={};
        
        parInfoL={};
        constrInfoL={};
        varInfoL={};
        
        sweepL={};
        sweepInfoL={};
        sweepParL={};
        plotSweep=true;
        
        report=              false;
        useReducedParList=   false;
    end
    
    properties (Dependent)
    end
    
    methods
        function obj=Proj(constP,basisP)
            if nargin > 0
                obj.basis=basisP;
                obj.const=constP;
                
            end
            obj.fitL={};
        end
        
        function obj=addFit(obj,fit)
            
            nFit=size(obj.fitL,1);
            obj.fitL{nFit+1,1}=fit;
        end
        
        %
        function L=exptL(obj)
            nFit=size(obj.fitL,1);
            L={};
            for iFit=1:nFit
                fit=obj.fitL{iFit};
                expt=fit.expt;
                newExpt=true;
                nExpt=size(L,1);
                for iExpt=1:nExpt
                    if isequal(expt,L{iExpt})
                        newExpt=false;
                        break;
                    end
                end
                if newExpt
                    L{nExpt+1,1}=expt;
                end
            end
        end
        
        %
        function simL=simL(proj)
            nFit=size(proj.fitL,1);
            simL={};
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                sim=fit.sim;
                newSim=true;
                nSim=size(simL,1);
                for iSim=1:nSim
                    if isequal(sim,simL{iSim})
                        newSim=false;
                        break;
                    end
                end
                if newSim
                    simL{nSim+1,1}=sim;
                end
            end
        end
        
        %
        function L=phaseL(obj)
            nFit=size(obj.fitL,1);
            L={};
            for iFit=1:nFit
                fit=obj.fitL{iFit};
                region=fit.region;
                nDomain=size(region.domainL,1);
                for iDomain=1:nDomain
                    nPhase=size(L,1);
                    newPhase=true;
                    phase=region.domainL{iDomain}.phase;
                    for iPhase=1:nPhase
                        if isequal(phase,L{iPhase})
                            newPhase=false;
                            break;
                        end
                    end
                    if newPhase
                        L{nPhase+1,1}=phase;
                    end
                end
            end
        end
        
        %
        function L=regionL(proj)
            nFit=size(proj.fitL,1);
            L={};
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                region=fit.region;
                nRegion=size(L,1);
                newRegion=true;
                for iRegion=1:nRegion
                    if isequal(region,L{iRegion})
                        newRegion=false;
                        break;
                    end
                end
                if newRegion
                    L{nRegion+1,1}=region;
                end
            end
            
        end
        
        %
        function L=domainL(proj)
            nFit=size(proj.fitL,1);
            L={};
            for iFit=1:nFit
                region=proj.regionL{iFit};
                nDomainP=size(region.domainL,1);
                for iDomain=1:nDomainP
                    nDomain=size(L,1);
                    newDomain=true;
                    domain=region.domainL{iDomain};
                    for iDomain=1:nDomain
                        if isequal(domain,L{iDomain})
                            newDomain=false;
                            break;
                        end
                    end
                    if newDomain
                        L{nDomain+1,1}=domain;
                    end
                end
            end
        end
        
        %
        function simL_usingRegion=simL_usingRegion(proj,region)
            simL_usingRegion={};
            nSim_usingRegion=0;
            
            simL=proj.simL();
            nSim=size(proj.simL,1);
            for iSim=1:nSim
                sim=simL{iSim};
                if isequal(region,sim.region)
                    simL_usingRegion{nSim_usingRegion+1,1}=sim;
                end
            end
        end
        
        function simL_usingPhase=simL_usingPhase(proj,phase)
            simL_usingPhase={};
            simL=proj.simL();
            nSim=size(proj.simL,1);
            for iSim=1:nSim
                phaseL=simL{iSim}.region.phaseL();
                nPhase=size(phaseL,1);
                for iPhase=1:nPhase
                    if isequal(phase,phaseL{iPhase})
                        nSim_usingPhase=size(simL_usingPhase,1);
                        simL_usingPhase{nSim_usingPhase+1,1}=simL{iSim};
                    end
                end
            end
        end
        
        function simL_usingDomain=simL_usingDomain(proj,domain)
            simL_usingDomain={};
            simL=proj.simL();
            nSim=size(proj.simL,1);
            for iSim=1:nSim
                domainL=simL{iSim}.region.domainL;
                nDomain=size(domainL,1);
                for iDomain=1:nDomain
                    if isequal(domain,domainL{iDomain})
                        nSim_usingDomain=size(simL_usingDomain,1);
                        simL_usingDomain{nSim_usingDomain+1,1}=simL{iSim};
                    end
                end
            end
        end
        
        %
        function [foundRegion region]=findRegion(proj,regionName)
            foundRegion=false;
            region=Region;
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                region=fit.region;
                if strcmp(region.name,regionName)
                    foundRegion=true;
                    return
                end
            end
        end
        
        %
        function [foundPhase phase]=findPhase(proj,phaseName)
            foundPhase=false;
            phase=Phase;
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                region=fit.region;
                nDomain=size(region.domainL,1);
                for iDomain=1:nDomain
                    phase=region.domainL{iDomain}.phase;
                    if strcmp(phase.name,phaseName)
                        foundPhase=true;
                        return;
                    end
                end
            end
        end
        
        %
        function [foundDomain domain]=findDomain(proj,regionName,phaseName)
            foundDomain=false;
            domain=Domain;
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                region=fit.region;
                if region.name==regionName
                    nDomain=size(region.domainL,1);
                    for iDomain=1:nDomain
                        domain=region.domainL{iDomain};
                        phase=domain.phase;
                        if strcmp(phase.name,phaseName)
                            foundDomain=true;
                            return;
                        end
                    end
                end
            end
        end
        
        %
        function ordL=ordL(proj)
            ordL={};
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                region=fit.region;
                nDomainP=size(region.domainL,1);
                for iDomain=1:nDomainP
                    domain=region.domainL{iDomain};
                    phase=domain.phase;
                    nOrd=size(ordL,1);
                    if strcmp(phase.type,'ordered')
                        ordL{nOrd+1,1}=phase;
                    end
                end
            end
        end
        
        %
        function proj=orientRegions(proj)
            const=proj.const;
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                region=fit.region;
                diff=fit.diff;
                region.orient(const,diff);
            end
        end
        
        %
        function proj=getSimBeamL(proj)
            
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit}.getSimBeamL(proj.const);
            end
        end
        
        %
        function proj=readExptBeamL(proj)
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                proj.fitL{iFit}.readExptBeamL();
            end
        end
        
        function proj=createFitLists(proj)
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                proj.fitL{iFit}.createBeamLists();
            end
        end
        
        %Create plot lists
        function proj=createPlotLists(proj)
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                fit.expt.createPlotList(fit.luL(:,1),fit.plotFitBeamsOnly,fit.diff.gMax);
                fit.sim.createPlotList(fit.luL(:,2),fit.plotFitBeamsOnly,fit.diff.gMax);
            end
        end
        
        %
        function parValV=getParValV(proj)
            nPar=size(proj.parInfoL,1);
            parValV=zeros(nPar,1);
            for iPar=1:nPar
                info=proj.parInfoL{iPar};
                switch info.type
                    case 1
                        x=log(info.region.thick);
                        parValV(iPar,1)=x;
                        
                    case 2
                        x=info.phase.eta/info.phase.etaMax;
                        parValV(iPar,1)=-log(1/x-1);
                        
                    case 3
                        x=log(info.domain.relVol);
                        parValV(iPar,1)=x;
                end
            end
        end
        
        %
        function proj=setParValV(proj,parValV)
            nPar=size(parValV,1);
            for iPar=1:nPar
                info=proj.parInfoL{iPar};
                switch info.type
                    case 1  %thickness
                        y=exp(parValV(iPar));
                        %yPrev=info.region.thick;
                        info.region.thick=y;
                        
                    case 2% eta
                        y=1/(1+exp(-parValV(iPar)));
                        parNew=info.phase.etaMax*y;
                        parPrev=info.phase.eta;
                        if parNew~=parPrev
                            nSim=size(info.simL,1);
                            for iSim=1:nSim
                                info.simL{iSim}.needsDiag=true;
                            end
                            info.phase.eta=parNew;
                        end
                        
                    case 3  %   relative volume
                        y=exp(parValV(iPar));
                        parNew=y;
                        parPrev=info.domain.relVol;
                        if parNew~=parPrev
                            nSim=size(info.simL,1);
                            for iSim=1:nSim
                                info.simL{iSim}.needsDiag=true;
                            end
                            info.domain.relVol=parNew;
                        end
                end
            end
            
            proj.adjust();
        end
        
        function parTolV=getParTolV(proj)
            nPar=size(proj.parInfoL,1);
            parTolV=zeros(nPar,1);
            
            for iPar=1:nPar
                info=proj.parInfoL{iPar};
                switch info.type
                    case 1
                        y=info.region.thick+1e-8;
                        parTolV(iPar)=1e-8/y;
                        
                    case 2
                        y=info.phase.eta*(info.phase.etaMax-info.phase.eta)+1e-8;
                        parTolV(iPar)=1e-8*info.phase.etaMax/y;
                        
                    case 3
                        y=info.domain.relVol+1e-8;
                        parTolV(iPar)=1e-8/y;
                end
            end
            
        end
        
        %
        function getParInfoL(proj)
            nPar=size(proj.parL,1);
            proj.parInfoL={};
            
            jPar=0;
            for iPar=1:nPar
                parDesc=proj.parL{iPar};
                info={};
                foundPar=false;
                switch parDesc{1}
                    case 'thickness'
                        info.type=1;
                        regionName=parDesc{2};
                        [res info.region]=proj.findRegion(regionName);
                        %[res info.level]=parDesc{3};%'main' or 'sub'
                        if res
                            foundPar=true;
                        else
                            "Region "+regionName+" for parameter "+parDesc{1}+" not found."
                        end
                        
                    case 'eta'
                        if ~proj.useReducedParList
                            info.type=2;
                            phaseName=parDesc{2};
                            [res info.phase]=proj.findPhase(phaseName);
                            if res
                                info.simL=proj.simL_usingPhase(info.phase);
                                foundPar=true;
                            else
                                "Phase "+phaseName+" for parameter "+parDesc{1}+" not found."
                            end
                        end
                        
                    case 'volume'
                        if ~proj.useReducedParList
                            info.type=3;
                            regionName=parDesc{2};
                            phaseName=parDesc{3};
                            [res info.domain]=proj.findDomain(regionName,phaseName);
                            if res
                                info.simL=proj.simL_usingDomain(info.domain);
                                foundPar=true;
                            else
                                "Domain in region "+regionName+" with phase "+phaseName+" for parameter "+parDesc{1}+" not found."
                            end
                            
                        end
                        
                    otherwise
                        "Parameter "+parDesc{1}+" not found."
                end
                if foundPar
                    jPar=jPar+1;
                    proj.parInfoL{jPar,1}=info;
                end
                
            end
            
        end
        
        %
        function proj=getVarInfoL(proj)
            nVar=size(proj.varL,1);
            proj.varInfoL={};
            
            jVar=0;
            for iVar=1:nVar
                varDesc=proj.varL{iVar,1};
                foundVar=false;
                info={};
                switch varDesc{1,1}
                    case 'thickness'
                        info.type=1;
                        regionName=varDesc{2,1};
                        [foundRegion info.region]=proj.findRegion(regionName);
                        if foundRegion
                            foundVar=true;
                        end
                        
                    case 'eta'
                        info.type=2;
                        phaseName=varDesc{2,1};
                        [foundPhase info.phase]=proj.findPhase(phaseName);
                        if foundPhase
                            foundVar=true;
                        end
                        
                    case 'volume fraction'
                        info.type=3;
                        regionName=varDesc{2,1};
                        phaseName=varDesc{3};
                        [foundDomain info.domain]=proj.findDomain(regionName,phaseName);
                        if foundDomain
                            foundVar=true;
                        end
                        
                    case 'eta net'
                        info.type=4;
                        info.ordL=proj.ordL();
                        foundVar=true;
                        
                    otherwise
                        "Variable "+varDesc{1,1}+" not found."
                        
                end
                if foundVar
                    jVar=jVar+1;
                    proj.varInfoL{jVar,1}=info;
                end
            end
            
        end
        
        %
        function proj=getConstrInfoL(proj)
            nConstr=size(proj.constrL,1);
            proj.constrInfoL={};
            
            jConstr=0;
            for iConstr=1:nConstr
                constrDesc=proj.constrL{iConstr};
                info={};
                foundVar=false;
                switch constrDesc{1}
                    case 'eta equal'
                        info.type=1;
                        phase1Name=constrDesc{2};
                        phase2Name=constrDesc{3};
                        [foundPhase1 info.phase1]=proj.findPhase(phase1Name);
                        [foundPhase2 info.phase2]=proj.findPhase(phase2Name);
                        if foundPhase1&&foundPhase2
                            if ~isequal(info.phase1,info.phase2)
                                foundVar=true;
                            end
                        end
                        
                    case 'eta/volume ratio constant'
                        info.type=2;
                        info.ordL=proj.makeOrdL();
                        foundVar=true;
                        
                    case 'eta max'
                        info.type=3;
                        phaseName=constrDesc{2};
                        info.phase=Phase;
                        [foundPhase info.phase]=proj.findPhase(phaseName)
                        if foundPhase
                            foundVar=true;
                        end
                        
                    otherwise
                        "Constraint "+constrDesc{1}+" not found."
                end
                if foundVar
                    jConstr=jConstr+1;
                    proj.constrInfoL{jConstr,1}=info;
                end
            end
            
        end
        
        %
        function getSweepInfoL(proj)
            nSweep=size(proj.sweepL,1);
            proj.sweepInfoL={};
            
            nSweep=size(proj.sweepL,1);
            jSweep=0;
            for iSweep=1:nSweep
                sweepDesc=proj.sweepL{iSweep};
                info={};
                foundSweep=false;
                switch sweepDesc{1}
                    case 'thickness'
                        info.type=1;
                        regionName=sweepDesc{2};
                        [res info.region]=proj.findRegion(regionName);
                        %[res info.level]=parDesc{3};%'main' or 'sub'
                        if res
                            info.rangeL=sweepDesc{3};
                            foundSweep=true;
                        else
                            "Region "+regionName+" for sweep "+sweepDesc{1}+" not found."
                        end
                        
                    case 'eta'
                        info.type=2;
                        phaseName=sweepDesc{2};
                        [res info.phase]=proj.findPhase(phaseName);
                        if res
                            info.rangeL=sweepDesc{3};
                            info.simL=proj.simL_usingPhase(info.phase);
                            foundSweep=true;
                        else
                            "Phase "+phaseName+" for sweep "+sweepDesc{1}+" not found."
                        end
                        
                    case 'volume'
                        info.type=3;
                        regionName=sweepDesc{2};
                        phaseName=sweepDesc{3};
                        [res info.domain]=proj.findDomain(regionName,phaseName);
                        if res
                            info.rangeL=sweepDesc{4};
                            info.simL=proj.simL_usingDomain(info.domain);
                            foundSweep=true;
                        else
                            "Domain in region "+regionName+" with phase "+phaseName+" for sweep "+sweepDesc{1}+" not found."
                        end
                        
                    otherwise
                        "Sweep "+sweepDesc{1}+" not found."
                end
                if foundSweep
                    jSweep=jSweep+1;
                    proj.sweepInfoL{jSweep,1}=info;
                end
                
            end
            
        end
        
        %setSweepVal(iSweepL);
        function proj=setSweepVal(proj,iSweepL)
            
            nSweep=size(iSweepL,1);
            for iSweep=1:nSweep
                info=proj.sweepInfoL{iSweep};
                switch info.type
                    case 1  %thickness
                        valPrev=info.region.thick;
                        valNew=info.rangeL(iSweepL(iSweep));
                        info.region.thick=valNew;
                        
                    case 2 %eta
                        valPrev=info.phase.eta;
                        valNew=info.rangeL(iSweepL(iSweep));
                        if valNew~=valPrev
                            nSim=size(info.simL,1);
                            for iSim=1:nSim
                                info.simL{iSim}.needsDiag=true;
                            end
                            info.phase.eta=valNew;
                        end
                        
                    case 3  %relative volume
                        valPrev=info.domain.relVol;
                        valNew=info.rangeL(iSweepL(iSweep));
                        if valNew~=valPrev
                            nSim=size(info.simL,1);
                            for iSim=1:nSim
                                info.simL{iSim}.needsDiag=true;
                            end
                            info.domain.relVol=valNew;
                        end
                end
            end
            
            proj.adjust();
        end
        
        %
        function proj=adjust(proj)
            
            %Apply constraints
            nConstr=size(proj.constrInfoL,1);
            for iConstr=1:nConstr
                info=proj.constrInfoL{iConstr};
                switch info.type
                    case 1% same order parameter
                        info.phase1.eta=info.phase2.eta;
                        
                    case 2%Order parameter proportional to phase volume
                        nOrd=size(info.ordL,1);
                        %Find order-parameter and relative-volume totals
                        etaTot=0;relVolTot=0;
                        for iOrd=1:nOrd
                            domain=info.ordL{iOrd};
                            etaTot=etaTot+domain.phase.eta;
                            relVolTot=relVolTot+domain.relVol;
                        end
                        
                        fac=etaTot/relVolTot;
                        %Scale order parameters
                        for iOrd=1:nOrd
                            domain=info.ordL{iOrd};
                            domain.phase.eta=fac*domain.relVol;
                        end
                        
                    case 3% max order parameter
                        info.phase.eta=info.phase.etaMax;
                end
            end
            
            %Adjust phases for order parameter.
            nPhase=size(proj.phaseL(),1);
            for iPhase=1:nPhase
                phase=proj.phaseL{iPhase};
                phase.adjust();
                
            end
            
            %Adjust regions for wedge and coherence.
            nRegion=size(proj.regionL(),1);
            for iRegion=1:nRegion
                proj.regionL{iRegion}.adjust();
            end
            
        end
        
        %Find error vector
        function X_V=getX_V(proj)
            X_V=[];
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                fit=proj.fitL{iFit};
                X_V=[X_V;fit.getX_V(proj.const)];
            end
        end
        
        %
        function proj=plotExpts(proj)
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                proj.fitL{iFit}.plotExpt(proj.const);
            end
        end
        
        %
        function proj=plotSims(proj)
            nFit=size(proj.fitL,1);
            for iFit=1:nFit
                proj.fitL{iFit}.plotSim(proj.const);
            end
        end
        
        %
        function proj=getSimIntensL(proj)
            nFit=size(proj.fitL,1);
            intensL=[];
            for iFit=1:nFit
                intensL=[intensL proj.fitL{iFit}.getIntensL(proj.const)];
            end
        end
        
        function varValV=getVarValV(proj)
            
            nVar=size(proj.varInfoL,1);
            varValV=zeros(nVar,1);
            
            nRegion=size(proj.regionL(),1);
            for iRegion=1:nRegion
                proj.regionL{iRegion}.getVolFracV();
            end
            
            for iVar=1:nVar
                info=proj.varInfoL{iVar};
                switch info.type
                    case 1%region: thickness
                        varValV(iVar,1)=info.region.thick;
                        
                    case 2%phase: eta
                        varValV(iVar,1)=info.phase.eta;
                        
                    case 3%domain: volume fraction
                        volFrac=info.domain.relVol/info.domain.region.relVolTot;
                        varValV(iVar,1)=volFrac;
                        
                    case 4%proj: eta net
                        eta=0;
                        nRegion=size(proj.regionL(),1);
                        for iRegion=1:nRegion
                            eta=eta+proj.regionL{iRegion,1}.getEta();
                        end
                        varValV(iVar,1)=eta/nRegion;
                end
            end
        end
        
        %
        function delVarV=getDelVarV(proj,parValV,dParV,delParV)
            
            nVar=size(proj.varInfoL,1);
            nPar=size(proj.parInfoL,1);
            JM=zeros(nVar,nPar);
            for iPar=1:nPar
                parP_ValV=parValV;
                parP_ValV(iPar)=parP_ValV(iPar)+dParV(iPar)/2;
                varValP_V=proj.setParValV(parP_ValV).getVarValV();
                
                parM_ValV=parValV;
                parM_ValV(iPar)=parM_ValV(iPar)-dParV(iPar)/2;
                varValM_V=proj.setParValV(parM_ValV).getVarValV();
                
                for iVar=1:nVar
                    JM(iVar,iPar)=(varValP_V(iVar)-varValM_V(iVar))/dParV(iPar);
                end
            end
            delVarV=zeros(nVar,1);
            for iVar=1:nVar
                for iPar=1:nPar
                    delVarV(iVar)=(delVarV(iVar)^2+(JM(iVar,iPar)*delParV(iPar))^2)^0.5;
                end
            end
            
        end
        
        %
        function J_M=getJacobiM(proj,parValV,dParV)
            nPar=size(parValV,1);
            J_M=zeros(0,nPar);
            for iPar=1:nPar
                parP_ValV=parValV;
                parP_ValV(iPar)=parP_ValV(iPar)+dParV(iPar)/2;
                proj.setParValV(parP_ValV);
                XpV=proj.getX_V();
                
                parM_ValV=parValV;
                parM_ValV(iPar)=parM_ValV(iPar)-dParV(iPar)/2;
                proj.setParValV(parM_ValV);
                %FOMmV=get_FOM_V_thickness(proj);
                XmV=proj.getX_V();
                
                nBeamTot=size(XmV,1);
                for iBeam=1:nBeamTot
                    J_M(iBeam,iPar)=(XpV(iBeam)-XmV(iBeam))/dParV(iPar);
                end
            end
            
        end
        
        %
        function proj=reportVals0(proj,sLabel,redX2,varValV)
            sRes=sLabel+compose(": reduced X^2 = %10.1f",redX2);sFormat="%s";
            nVar=size(varValV,1);
            for iVar=1:nVar
                resV=[varValV(iVar)];
                info=proj.varInfoL{iVar};
                sResP="";
                switch info.type
                    case 1
                        sResP=info.region.name;
                        sResP=sResP+compose(": thick (nm) = %5.3f",resV);
                        
                    case 2
                        sResP=info.phase.name;
                        sResP=sResP+compose(": eta = %8.6f",resV);
                        
                    case 3%'volume fraction'
                        sResP=info.domain.region.name+", "+info.domain.phase.name;
                        sResP=sResP+compose(": volFrac = %8.6f",resV);
                        
                    case 4%'eta net'
                        sResP=compose("net: eta = %8.6f",resV);
                end
                sRes=[sRes sResP];
                sFormat=sFormat+"\n%s";
            end
            sprintf(sFormat,sRes)
        end
        
        %
        function proj=reportVals(proj,iIter,redX2,varValV,delVarV)
            sRes=compose("iteration %3.0f: reduced X^2 = %10.5f",[iIter redX2]);sFormat="%s";
            nVar=size(varValV,1);
            for iVar=1:nVar
                resV=[varValV(iVar) delVarV(iVar)];
                info=proj.varInfoL{iVar};
                sResP="";
                switch info.type
                    case 1
                        sResP=info.region.name;
                        sResP=sResP+compose(": thick (nm) = %5.3f(%5.3f)",resV);
                        
                    case 2
                        sResP=info.phase.name;
                        sResP=sResP+compose(": eta = %8.6f(%8.6f)",resV);
                        
                    case 3%'volume fraction'
                        sResP=info.domain.region.name+", "+info.domain.phase.name;
                        sResP=sResP+compose(": volFrac = %8.6f(%8.6f)",resV);
                        
                    case 4%'eta net'
                        sResP=compose("net: eta = %8.6f(%8.6f)",resV);
                end
                sRes=[sRes sResP];
                sFormat=sFormat+"\n%s";
            end
            sprintf(sFormat,sRes)
        end
        
        %
        function proj=reportSweep(proj,sDesc,redX2,iSweepL)
            sRes=sDesc+compose(", reduced X^2 = %10.5f",redX2);sFormat="%s";
            nSweep=size(iSweepL,1);
            for iSweep=1:nSweep
                info=proj.sweepInfoL{iSweep};
                val=info.rangeL(iSweepL(iSweep));
                sResP="";
                switch info.type
                    case 1
                        sResP=info.region.name;
                        sResP=sResP+compose(": thick (nm) = %5.3f",val);
                        
                    case 2
                        sResP=info.phase.name;
                        sResP=sResP+compose(": eta = %8.6f",val);
                        
                    case 3%'volume fraction'
                        sResP=info.domain.region.name+", "+info.domain.phase.name;
                        sResP=sResP+compose(": volume = %8.6f",val);
                end
                sRes=[sRes sResP];
                sFormat=sFormat+"\n%s";
            end
            sprintf(sFormat,sRes)
        end
        
        %Prepare sims for simulation
        function proj=prepare(proj)
            
            phaseL=proj.phaseL();
            nPhase=size(phaseL,1);
            ordL=proj.ordL();
            nOrd=size(ordL,1);
            etaL=zeros(nOrd,1);
            for iOrd=1:nOrd
                ord=ordL{iOrd};
                etaL(iOrd)=ord.eta;
                ord.eta=ord.etaMax;
            end
            
            proj.adjust();
            
            simL=proj.simL();
            nSim=size(simL,1);
            for iSim=1:nSim
                simL{iSim}.prepare(proj.const);
            end
            
            %Now restore order params
            for iOrd=1:nOrd
                ord=ordL{iOrd};
                ord.eta=etaL(iOrd);
            end
            proj.adjust();
            
        end
        
        %Summary of fit values
        function proj=summarize(proj)
            
            simL=proj.simL();
            nSim=size(simL,1);
            for iSim=1:nSim
                simL{iSim}.summarize();
            end
            
        end
        
        % Orient bases, read peak files, create simulations
        function proj=init(proj)
            
            proj.orientRegions();
            proj.readExptBeamL();
            proj.getSimBeamL();
            proj.createFitLists();
            
            proj.createPlotLists();
            proj.plotExpts();
            
            proj.getConstrInfoL();
            proj.getVarInfoL();
            
            proj.prepare();
            proj.summarize();
        end
        
        % General fit
        function res=fit(proj)
            proj.getParInfoL();
            nPar=size(proj.parInfoL,1);
            parValV=proj.getParValV();
            parTolV=zeros(nPar,1);
            
            varValV=proj.getVarValV();
            X_bestV=proj.getX_V();
            X2_best=norm(X_bestV)^2;
            
            nData=size(X_bestV,1);
            nPar=size(proj.parL,1);
            X2_redBest=X2_best/(nData-nPar);
            
            if proj.report
                proj.reportVals0("initial",X2_redBest,varValV);
            end
            
            proj.plotSims();
            
            iIter=1;
            nu=0.1;
            done=(nPar==0);
            newBest=true;
            %needJM=true;
            while (iIter<=proj.nIters)&&(~done)% iteration loop
                parTolV=proj.getParTolV();
                dParV=parTolV*1e3;
                
                %
                if newBest
                    J_M=proj.getJacobiM(parValV,dParV);
                end
                
                %
                got1=false;
                iJ_M=getDampedPseudoInverseLeft(J_M,nu);
                par1IncV=(-(iJ_M*X_bestV));
                %adjustParIncV(par1IncV)
                par1ValV=parValV+par1IncV;
                proj.setParValV(par1ValV);
                var1ValV=proj.getVarValV();
                X_1V=proj.getX_V();
                X2_1=norm(X_1V)^2;
                %proj.reportVals0("nu",FOM1,var1ValV);
                got1=true;
                
                %
                got2=false;
                iJ_M=getDampedPseudoInverseLeft(J_M,nu/2);
                par2IncV=(-(iJ_M*X_bestV));
                par2ValV=parValV+par2IncV;
                proj.setParValV(par2ValV);
                var2ValV=proj.getVarValV();
                
                X_2V=proj.getX_V();
                X2_2=norm(X_2V)^2;
                %proj.reportVals0("nu/2",FOM2,var2ValV);
                got2=true;
                
                %
                if got1&&got2
                    if (X2_1<X2_best)||(X2_2<X2_best)
                        if X2_1<X2_2
                            parValV=par1ValV;
                            dParV=par1IncV/1e3;
                            X2_best=X2_1;
                            X_bestV=X_1V;
                            
                        else
                            parValV=par2ValV;
                            dParV=par2IncV/1e3;
                            X2_best=X2_2;
                            Xbest_V=X_2V;
                            nu=nu/2;
                        end
                        newBest=true;
                        %needJM=true;
                    else
                        nu=nu*10;
                        newBest=false;
                        if nu>1e6
                            compose("Unable to improve fit.")
                            done=true;
                        end
                        %needJM=false;
                    end
                else
                    done=true;
                end
                
                %dParV(dParV<parTolV)=parTolV(dParV<parTolV);
                if newBest
                    proj.setParValV(parValV);
                    varValV=proj.getVarValV();
                    
                    %Evaluate parameter errors
                    nData=size(X_bestV,1);
                    nPar=size(proj.parL,1);
                    
                    iCovarM=getCovarianceInverseM(J_M);
                    delParV=sqrt(diag(iCovarM));
                    
                    %Find variable errors
                    delVarV=proj.getDelVarV(parValV,dParV,delParV);
                    
                    %Report iteration
                    X2_redBest=X2_best/(nData-nPar);
                    if proj.report
                        proj.reportVals(iIter,X2_redBest,varValV,delVarV);
                    end
                    
                    proj.plotSims();
                    iIter=iIter+1;
                end%if newBest
            end%while (iIter<=nIters)&&(~done)
            res.X2_best=X2_best;
            res.X2_redBest=X2_redBest;
            res.parValV=parValV;
            res.varValV=varValV;
           
        end
        %
        function proj=fastFit(proj)
            
            for iIter=1:10
                compose("fit")
                proj.nIters=3;
                proj.useReducedParList=false;
                proj.fit();
                
                compose("reduced fit")
                proj.nIters=10;
                proj.useReducedParList=true;
                proj.fit();
            end
        end
        
        %
        function res=sweep(proj)
            proj.report=false;
            
            fullParL=proj.parL;
            proj.parL=proj.sweepParL;
            proj.getParInfoL();
            parInitValV=proj.getParValV();
           
            proj.getSweepInfoL();
            nSweep=size(proj.sweepInfoL,1);
            iSweepL=ones(nSweep,1);
            iSweepBestL=ones(nSweep,1);
            
            sizeM=[];
            for iSweep=1:nSweep
                sizeM=[sizeM size(proj.sweepInfoL{iSweep}.rangeL,1)];
            end
            if(nSweep==1)
                X2_redM=zeros(sizeM,1);
            else
                X2_redM=zeros(sizeM);
            end
            nPointTot=1;
            for iSweep=1:nSweep
                nPointTot=nPointTot*sizeM(1,iSweep);
            end
            
            first=true;
            nPoint=0;
            done=false;
            while ~done
                proj.setSweepVal(iSweepL);
                proj.setParValV(parInitValV);
                fitRes=proj.fit();
               
                X2_M(nPoint+1)=fitRes.X2_best;
                X2_redM(nPoint+1)=fitRes.X2_redBest;
                X2_redCurr=fitRes.X2_redBest;
                X2_curr=fitRes.X2_best;
                 
                %Keep track of best FOM and values
                if first
                    X2_redBest=X2_redCurr;
                    X2_best=X2_curr;
                    iSweepBestL=iSweepL;
                end
                if X2_redCurr<X2_redBest
                    X2_redBest=X2_redCurr;
                    X2_best=X2_curr;
                    iSweepBestL=iSweepL;
                end
                first=false;
                
                nPoint=nPoint+1;
                fProg=100*nPoint/nPointTot;
                %compose("----progress = %5.2f",fProg)
                proj.reportSweep("current",X2_redCurr,iSweepL);
                 
                
                for iSweep=1:nSweep
                    iSweepL(iSweep)=iSweepL(iSweep)+1;
                    iSweepMax=size(proj.sweepInfoL{iSweep}.rangeL,1);
                    if iSweepL(iSweep)>iSweepMax
                        iSweepL(iSweep)=1;
                        if iSweep==nSweep
                            done=true;
                        end
                    else
                        break
                    end
                end
                
                if iSweep>1
                    compose("----progress = %5.2f",fProg)
                    proj.reportSweep("best",X2_redBest,iSweepBestL);
                end
            end
            
            
            res.X2_best=X2_best;
            res.X2_redBest=X2_redBest;
            res.X2_redM=X2_redM;
            res.X2_M=X2_M;
            res.iSweepBestL=iSweepBestL;
            %proj.sweepInfoL
            if proj.plotSweep&&(nSweep==1)
                res.xL=proj.sweepInfoL{1}.rangeL;%rows
                figure
                plot(res.xL,res.X2_M)
                set(gca,'xscale','log')
            end
            
            if proj.plotSweep&&(nSweep>1)
                res.yL=proj.sweepInfoL{1}.rangeL;%columns
                res.xL=proj.sweepInfoL{2}.rangeL;%rows
                figure
                surf(res.xL,res.yL,res.X2_redM)
                set(gca,'yscale','log')
            end
            
            proj.parL=fullParL;
        end
        
        %Performs a sweep to find best global fit using sweepParL.
        %Then fits usint parL.
        function sweepFit(proj)
            
            proj.report=false;
            res=proj.sweep();
            
            proj.setSweepVal(res.iSweepBestL);
            proj.report=true;
            proj.fastFit();
        end
        
        
    end
end

