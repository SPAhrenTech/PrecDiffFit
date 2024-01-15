classdef a_nm
    properties
        
        aAlP=                       0.54510;%nm
        aGaP=                       0.54512;%nm
        aInP=                       0.58686;%nm

        aAlAs=                      0.54510;%nm
        aGaAs=                      0.54512;%nm
        aInAs=                       0.58686;%nm

        aAlSb=                      0.54510;%nm
        aGaSb=                      0.54512;%nm
        aInSb=                      0.58686;%nm
    end
    methods
        function a = get_a(obj,w,x,y,z)
            xAl=w*x;
            xGa=(1-w)*x;
            xIn=1-x;

            yP=1-y;
            yAs=y*(1-z);
            ySb=y*z;

            aAl=yP*obj.aAlP+yAs*obj.aAlAs+ySb*obj.aAlSb;
            aGa=yP*obj.aGaP+yAs*obj.aGaAs+ySb*obj.aGaSb;
            aIn=yP*obj.aInP+yAs*obj.aInAs+ySb*obj.aInSb;

            a=xAl*aAl+xGa*aGa+xIn*aIn;
        end
    end
end
