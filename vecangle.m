function result = vecangle(x1,y1,z1,x2,y2,z2)
if nargin<3
    z2=y1(3);
    z1=x1(3);
    y2=y1(2);
    x2=y1(1);
    y1=x1(2);
    x1=x1(1);
end
result=acosd((x1.*x2+y1.*y2+z1.*z2)./(sqrt(x1.^2+y1.^2+z1.^2).*sqrt(x2.^2+y2.^2+z2.^2)));