classdef Diff<handle
    properties
        EKeV=                 200;%beam energy <- Input
        lambda;               
        T_K=                  300;%Not currently used.

        phi=                  40;%beam tilt (mrad)
        nSteps=               180;%# of precession steps
        N_ZoneL=               [0];%zone index <- Input
        sgRange=              0.02;%Buerger precession filter range (1/nm)
        mode=                 'max';%maximum: 'max',integrated: 'int'
        doBuerger=            false;%Use Buerger mask for precession
        nFoilV;
        kV;                  %incident beam wave vector (1/nm)

        %
        tExposure=         4;%
        fContrast=         0;%
        camlen=            10000;%pix
        pixX=              512;%image size X
        pixY=              512;%image size Y
        const;
        
        gMax=               1e3;%(1/nm)
        autoExp=            false;
        
        doRealPlot=           false;
    end
    
    methods
        function diff=Diff(const,EKeV)
            if nargin>0
                diff.const=const;
                diff.EKeV=EKeV;
                diff.nFoilV=const.nZ;
                diff.lambda=const.hc/sqrt(EKeV*(EKeV+2*const.m0c2));
                diff.kV=const.nZ/diff.lambda;
            end
        end
     
        function img=plot(diff,const,region,intensL,H_refL,beamL)

            camConst=       diff.camlen*diff.lambda;

            x=              2*((1:diff.pixX)-diff.pixX/2)/camConst;
            y=              -2*((1:diff.pixY)-diff.pixY/2)/camConst;  

            gL=(region.B_refM'*H_refL')';

            [X,Y]=meshgrid(x,y);

            %image
            pixX=   size(X,2);
            pixY=   size(X,1);
            intensM=zeros(pixX,pixY);%whole pattern

            nBeams=size(beamL,1);
            intensMax=0;intensMin=1;
            for ii=1:nBeams
                if intensL(beamL(ii))>intensMax
                    intensMax=intensL(ii);
                end
                if intensL(beamL(ii))<intensMin
                   intensMin=intensL(beamL(ii));
                end
            end

            intens0=intensMax;
            A0=1;
            rad0=0.1;%/(vRef^(1/3))/4;
            for ii=1:nBeams
                gV=gL(beamL(ii),:);
                gX=dot(gV,const.nX);
                gY=dot(gV,const.nY);
                intens=(intensL(beamL(ii)));%+10^tExposure);
                f=(intens/intens0)^(1/3);
                A=(f^2)*A0;
                rad=(f^0.5)*rad0;
                intensM=intensM+A*exp(-((X-gX).^2+(Y-gY).^2)/rad^2);
              
            end
           img=(10^diff.tExposure)*intensM;    
           if diff.autoExp
               Itot=sum(p);
               img=255.*p./Itot;
           end
        end

      %
      function img=realPlot(diff,const,region,intensL,H_refL,beamL)

       camConst=       diff.camlen*diff.lambda;

        x=              2*((1:diff.pixX)-diff.pixX/2)/camConst;
        y=              -2*((1:diff.pixY)-diff.pixY/2)/camConst;  

        gL=(region.B_refM'*H_refL')';

        [X,Y]=meshgrid(x,y);

        %image
        pixX=   size(X,2);
        pixY=   size(X,1);
        intensM=zeros(pixX,pixY);%whole pattern

        nBeams=size(beamL,1);
        intensMax=0;intensMin=1;
        for ii=1:nBeams
            if intensL(beamL(ii))>intensMax
                intensMax=intensL(ii);
            end
            if intensL(beamL(ii))<intensMin
               intensMin=intensL(beamL(ii));
            end
        end

        intens0=intensMax;
        A0=1;
        rad0=0.25;%/(vRef^(1/3))/4;
        for ii=1:nBeams
            gV=gL(beamL(ii),:);
            gX=dot(gV,const.nX);
            gY=dot(gV,const.nY);
            intens=(intensL(beamL(ii)));%+10^tExposure);
            A=intens;
            rad=rad0;
            
            f=(intens/intens0)^(1/2);
            Ap=A/f;
            radp=(f^0.5)*rad;
            intensM=intensM+Ap*exp(-((X-gX).^2+(Y-gY).^2)/radp^2);

        end
       img=(10^diff.tExposure)*intensM;    
      end

    end

    
end

