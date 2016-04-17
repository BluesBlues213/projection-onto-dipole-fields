function [kdPole] = dipole2D(X)

xx=X(1);
yy=X(2);
kdPole=zeros([xx yy]);

for y=1:yy;
    for x=1:xx;
        r=(x-(xx+1)/2)^2+(y-(yy+1)/2)^2;
        kdPole(x,y)=1/3 -(x-(xx+1)/2)^2/r;
    end;
end;

kdPole((xx)/2,(yy)/2)=0;
