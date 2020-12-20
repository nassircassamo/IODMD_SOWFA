function [xx,yy,zz]=resamplegrid(x,y,z, decimate)

    a=0;
    for k=1:decimate:length(x)
        a=a+1;
        xx(a)=x(k);
    end

    a=0;
    for k=1:decimate:length(y)
        a=a+1;
        yy(a)=y(k);
    end

    a=0;
    for k=1:decimate:length(z)
        a=a+1;
        zz(a)=z(k);
    end
    
    xx=transpose(xx);
    yy=transpose(yy);
    zz=transpose(zz);

end
