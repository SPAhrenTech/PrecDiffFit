classdef OrdPhase<Phase
   properties
    eta=0.15;
    occDisordL=[];
    occOrdL=[];
	dispDisordL=[];
      dispOrdL=[];
      dispV=[];
      dispDist=0;
      etaMax;
   end
    
	methods
      function obj=OrdPhase(name,etaMaxP)
        obj@Phase(name);
         if nargin > 0
            obj.type='ordered';
            obj.etaMax=etaMaxP;               
        end
        end

        function obj=etaFromPar(obj,parVal)
             y=1/(1+exp(-parVal));
             obj.eta=obj.etaMax*y;
        end

        function parVal=parFromEta(obj)
            x=obj.eta/specimen.etaMax;
            parVal=-log(1/x-1);
        end
        
        %Adjust for order parameter
        %Correct occupation list
        function phase=adjust(phase)
            basis=phase.basis;
            phase.occL=(1-phase.eta)*phase.occDisordL...
                +(phase.eta)*phase.occOrdL;

            %Correct coordinates
            dispL=(1-phase.eta)*phase.dispDisordL...
                +(phase.eta)*phase.dispOrdL;

            %Displace anions
            nrOrd=basis.A_refM'*phase.dispV;
            nrOrd=nrOrd/norm(nrOrd);
            dispV=phase.dispDist*nrOrd;
            delDL=dispV*dispL';
            dL=basis.A_refM'*phase.coordL+delDL;
            phase.adjCoordL=basis.B_refM*dL;           
        end
        
        %
        function obj=AlInP_n111(obj,xIn,eta)
            w=1;x=1-xIn;y=0;z=0;
            a=a_nm().get_a(w,x,y,z);
            xAl=1-xIn;

            obj.A_primM=       	a*[	-1 2 1
                                -2 1 1
                                -1 1 2]/2;%Vectors are rows

            obj.elemL=            {'Al' 'In' 'In' 'Al' 'P' 'P'};

            obj.siteL=            [1 1 2 2 3 4]';
            obj.coordL=           [	0 0 0     
                                    -4 4 4
                                    1 1 1
                                    -3 5 5]'/4;%(-1 1 1)Right variant

            obj.occDisordL=       [xAl xIn xIn xAl 1 1]';
            obj.occOrdL=          [1/2+xAl 1/2-xAl 1/2+xIn 1/2-xIn 1 1]';

            obj.dispDisordL=      [0 0 0 0]';%Relative displacement for each disordered site
            obj.dispOrdL=         [0 0 -1 1]';%Relative displacement for each orderedsite
            obj.dispV=            [-1 1 1]';%displacement direction
            obj.dispDist=         0.0119746047519545;%nm, for fully ordered, from VFF
            obj.eta=        		eta;
            obj.absorp=             0.1;
        end
        
        %
       function obj=AlInP_1n11(obj,xIn,eta)
            w=1;x=1-xIn;y=0;z=0;
            a=a_nm().get_a(w,x,y,z);
            xAl=1-xIn;
            
            obj.A_primM=         	a*[ 1 -2 1
                                        2 -1 1
                                        1 -1 2]/2;%Vectors are rows

            obj.elemL=            {'Al' 'In' 'In' 'Al' 'P' 'P'};

            obj.siteL=            [1 1 2 2 3 4]';

            obj.coordL=             [   0 0 0  
                                        4 -4 4
                                        1 1 1
                                        5 -3 5]'/4;%(1 -1 1) Left variant

            obj.occDisordL=         [xAl xIn xIn xAl 1 1]';
            obj.occOrdL=          [1/2+xAl 1/2-xAl 1/2+xIn 1/2-xIn 1 1]';

            obj.dispDisordL=        [0 0 0 0]';%Relative displacement for each disordered site
            obj.dispOrdL=            [0 0 -1 1]';%Relative displacement for each orderedsite
            obj.dispV=              [1 -1 1]';%displacement direction
            obj.dispDist=       	0.0119746047519545;%nm, for fully ordered, from VFF
            obj.eta=        		eta;
            obj.absorp=             0.1;

       end
    end
end
