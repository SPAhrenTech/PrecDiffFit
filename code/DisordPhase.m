classdef DisordPhase<Phase
   properties
   end
   
   methods
        function obj=DisordPhase(name)
            obj@Phase(name);
            if nargin>0
                obj.type='disordered';
            end           
        end
        
        function adjust(phase)
            adjust@Phase(phase);
        end
        
      %
     function obj=Si(obj)
        a=0.54;
  
        obj.A_primM=        	a*[	0 1 1
                                    1 0 1
                                    1 1 0]/2;%Vectors are rows

       obj.elemL=               {   'Si' 'Si'};

       obj.siteL=               [   1 2]';
       obj.coordL=              [	0 0 0     %
                                    1 1 1]'/4;%wrt Aref

       obj.occL=                [	1 1]';
     end
       
        
        %
     function obj=AlInP(obj,xIn)
        w=1;x=1-xIn;y=0;z=0;
        a=a_nm().get_a(w,x,y,z);
        xAl=1-xIn;

        obj.A_primM=        	a*[	0 1 1
                                    1 0 1
                                    1 1 0]/2;%Vectors are rows

       obj.elemL=               {'Al' 'In' 'P'};

       obj.siteL=               [    1  1 2]';
       obj.coordL=              [	0 0 0     %
                                    1 1 1]'/4;%wrt Aref

       obj.occL=                [1-xIn xIn 1]';
     end
     
     
   end
end
