% --------------------------------------------------------------------------------
% Structure factors - P. Ahrenkiel (2020)
% --------------------------------------------------------------------------------
function Ug=getUg(const,sim,H_ref)
 
    region=sim.region;
    basis=region.basis;
    diff=sim.diff;
	gV=basis.B_refM'*H_ref;
	s=norm(gV)/2;
    Ug_V=0;%Hermitian and non-hermitian combined
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
                *(1+1i*phase.absorp);%complex           
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