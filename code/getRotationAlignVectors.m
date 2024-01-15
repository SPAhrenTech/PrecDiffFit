function RM=getRotationAlignVectors(V,Vp)

    u1=V/norm(V);
    u3=cross(u1,Vp);
    if(norm(u3)==0)
        u3=cross(u1,[0 1 0]');
    end
    if(norm(u3)==0)
        u3=cross([0 0 1]',u1);
    end
    u3=u3/norm(u3);
    u2=cross(u3,u1); 
    R1=[u1 u2 u3];
    
    up1=Vp/norm(Vp);
    up3=u3;
    up2=cross(up3,up1);
    up2=up2/norm(up2);
    R2=[up1 up2 up3];
    
    R1i=inv(R1);
    RM=(R2*R1i);
end

