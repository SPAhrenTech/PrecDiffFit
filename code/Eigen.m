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
            
           nStrong=size(eigen.strongL,1);
           if nStrong<2
            eigen.strongL=[1:nBeam]';
            eigen.pertL=[];
            eigen.negL=[];
           end
        end
        
      %Get the Bloch waves
      function eigen=calculate(eigen,const,sim,iXsiM)
          
            nBeam=size(sim.H_refL,1);             
            nStrong=size(eigen.strongL,1);
            nPert=size(eigen.pertL,1);
            
            if nPert>0
                [col,~]=meshgrid(1:nPert,1:nStrong);
                iXsiBetheM=iXsiM(eigen.strongL,eigen.pertL,1);
                isgBetheL=1./eigen.sgL(eigen.pertL);
                bethe1M=isgBetheL(col).*iXsiBetheM;
                bethe2M=bethe1M*iXsiBetheM';
            end
            
            %Find Bethe potentials and excitation errors 
            iXsiStrong_pM=iXsiM(eigen.strongL,eigen.strongL,1);
            sgStrong_pL=eigen.sgL(eigen.strongL,1);
            if nPert>0
            	iXsiStrong_pM=iXsiStrong_pM-1/2*bethe2M;
                sgStrong_pL=sgStrong_pL-1/4*diag(bethe2M);
            end
            A_strongM=iXsiStrong_pM/2;
            A_strongM(1:(nStrong+1):nStrong*nStrong)=sgStrong_pL;
            % A_strongM(sub2ind(size(A_strongM),1:nStrong,1:nStrong))=sgStrong_pL;
              
            %Diagonalize
            [C_strongM,G_M]=eig(A_strongM);
            nBloch=nStrong;
            
            %Now find Bloch waves coeffs
            eigen.C_M=zeros(nBeam,nBloch);
            eigen.C_M(eigen.strongL,:)=C_strongM(:,:);
            
            if nPert>0
                C_pertM=-1/2*bethe1M'*C_strongM;
            	eigen.C_M(eigen.pertL,1:nBloch)=C_pertM;   
            end
            
            %Make Bloch coeffs zero for neglected reflections
            eigen.C_M(eigen.negL,1:nBloch)=0; 
              
            %absorption 
            ApM=iXsiM(:,:,2)/2;
            GpM=(eigen.C_M')*ApM*eigen.C_M;
 
            eigen.gammaL=diag(G_M)+1i*diag(GpM);
        end
    end
end