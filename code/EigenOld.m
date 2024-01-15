classdef Eigen<handle
	properties
       	gammaL=[];
        C_M=[];
       % pertM;%Matrix specifying beam status at each precession angle
       % usePert=0;
        strongL=[];%included in diagonalization
        pertL=[];%estimated by perturbation
        negL=[];%assumed zero
        sgL=[];%excitation errors
    end
    
	methods
        function eigen=Eigen(nBeams)
      		if nargin > 0
                eigen.gammaL=zeros(nBeams,1);
                eigen.C_M=zeros(nBeams,nBeams);
            end
        end
         
        %
        function eigen=prepare(eigen,sim,iXsiM)            
            diff=sim.diff;
            region=sim.region;
  
            nBeam=size(sim.H_refL,1);
            eigen.sgL=zeros(nBeam,1);
            for iBeam=1:nBeam
                H_ref=sim.H_refL(iBeam,:)';
                gV=region.B_refM'*H_ref;
                eigen.sgL(iBeam)=getExcErr(gV,diff.kV,diff.nFoilV);
            end
            
            eigen.strongL=[];
            eigen.pertL=[];
            eigen.negL=[];
            if sim.usePert
                for iBeam=1:nBeam
                    if abs(eigen.sgL(iBeam))<sim.sgMax
                       eigen.strongL=[eigen.strongL;iBeam];

                    else
                       invw=norm(iXsiM(iBeam,1,1)/eigen.sgL(iBeam));
                       if invw>sim.invwMin
                           eigen.pertL=[eigen.pertL;iBeam];
                       else
                           eigen.negL=[eigen.negL;iBeam];
                       end
                    end
                end
            else
                eigen.strongL=[1:nBeam]';
            end
        end
        
      %Get the Bloch waves
      function eigen=calculate(eigen,const,sim,iXsiM)

            diff=sim.diff;
            region=sim.region;
            k=1/diff.lambda;
            
            nBeam=size(sim.H_refL,1);             
            nStrong=size(eigen.strongL,1);
            nPert=size(eigen.pertL,1);
            nNeg=size(eigen.negL,1);
            
            %Find Bethe potentials and excitation errors
            A_strongM=zeros(nStrong,nStrong);
            
            %%iBeamL=eigen.strongL;
            for iStrong=1:nStrong
                iBeam=eigen.strongL(iStrong);
                
                sgp=eigen.sgL(iBeam,1);%Start with actual excitation error of strong beams
                for jPert=1:nPert%Adjust excitation error for strong beams
                    jBeam=eigen.pertL(jPert);
                    sgp=sgp-abs(iXsiM(iBeam,jBeam,1))^2/4/eigen.sgL(jBeam);
                end
                A_strongM(iStrong,iStrong)=sgp;
                %diag(A_strongM)=sgpL;
                
               %jBeam=eigen.strongL;
                for jStrong=iStrong+1:nStrong
                   jBeam=eigen.strongL(jStrong);
 
                   iXsip=iXsiM(iBeam,jBeam,1);%Start with actual coupling between strong beams
                   %kBeamL=eigen.pertL;
                  % iXsipM=iXsipM-iXsiM(iBeam,kBeam,1)*iXsiM(kBeam,jBeam,1)/2/eigen.sgL(kBeam,1);
                   for kPert=1:nPert%Adjust coupling between strong beams
                    kBeam=eigen.pertL(kPert);
                    iXsip=iXsip-iXsiM(iBeam,kBeam,1)*iXsiM(kBeam,jBeam,1)/2/eigen.sgL(kBeam);
                   end
                  % A_strongM=iXsipM/2;
                    A_strongM(iStrong,jStrong)=iXsip/2;
                    A_strongM(jStrong,iStrong)=iXsip'/2;
                 
                end                    
            end
            
            %Diagonalize
            [C_strongM,G_M]=eig(A_strongM);
            eigen.gammaL=diag(G_M);
            nBloch=nStrong;
            
            %Now find Bloch waves coeffs
            eigen.C_M=zeros(nBeam,nBloch);         
            for iStrong=1:nStrong
                iBeam=eigen.strongL(iStrong);
                eigen.C_M(iBeam,:)=C_strongM(iStrong,:);             
            end
         
            %eigen.C_M(iBeamL,:)=C_strongM(:,:);
            
            %Find Bloch coeff's for perturbed reflections
            C_pertM=zeros(nPert,nBloch);
            for iPert=1:nPert
               iBeam=eigen.pertL(iPert,1);
               for jBloch=1:nBloch
                    C_pert=0;
                    for kStrong=1:nStrong
                        kBeam=eigen.strongL(kStrong);
                        C_pert=C_pert-iXsiM(iBeam,kBeam,1)*eigen.C_M(kBeam,jBloch)/2/eigen.sgL(iBeam);
                   end
                 eigen.C_M(iBeam,jBloch)=C_pert;
                 C_pertM(iPert,jBloch)=C_pert;
                end
            end
                    
            %Make Bloch coeffs zero for neglected reflections
            for iNeg=1:nNeg
                iBeam=eigen.negL(iNeg);
                 for jBloch=1:nBloch
                     eigen.C_M(iBeam,jBloch)=0;   
                end
            end
            
            %absorption          
            ApM=zeros(nBeam,nBeam);
            for iBeam=1:nBeam
                 ApM(iBeam,iBeam)=iXsiM(iBeam,iBeam,2)/2;             
                for jBeam=iBeam+1:nBeam
                  ApM(iBeam,jBeam)=iXsiM(iBeam,jBeam,2)/2;
                  ApM(jBeam,iBeam)=iXsiM(jBeam,iBeam,2)/2;
                end
            end
            GpM=(eigen.C_M')*ApM*eigen.C_M;
            gammapL=diag(GpM);

            for iBloch=1:nStrong
                eigen.gammaL(iBloch,1)=eigen.gammaL(iBloch,1)+1i*gammapL(iBloch,1);
            end
       end
    end
end