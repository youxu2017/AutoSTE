function Z_plotShipXY(xo,yo,course)
    L = 1.6;
    B = 0.9;
%     画图案的中心区域
%   plt1 =  plot(xo,yo,'ro');
    xo1=xo+0.5*L*cos(course);
    yo1=yo+0.5*L*sin(course);
    fd=0:0.04:1;
    x=xo1+fd*sqrt((0.5*B)^2+(0.25*L)^2)*cos(course+pi-atan(2*B/L));
    y=yo1+fd*sqrt((0.5*B)^2+(0.25*L)^2)*sin(course+pi-atan(2*B/L));
    plt2 = plot(x,y,'-b');
    hold on
    x=xo1+fd*sqrt((0.5*B)^2+(0.25*L)^2)*cos(course-(pi-atan(2*B/L)));
    y=yo1+fd*sqrt((0.5*B)^2+(0.25*L)^2)*sin(course-(pi-atan(2*B/L)));
    plt3 = plot(x,y,'-b');
    hold on
    xo2=xo+0.5*B*sin(course)-0.5*L*cos(course);
    yo2=yo-0.5*L*sin(course)-0.5*B*cos(course);
    x=xo2+fd*0.75*L*cos(course);
    y=yo2+fd*0.75*L*sin(course);
    plt4 = plot(x,y,'-b');
    hold on
    x=xo2+fd*B*cos(course-1.5*pi);
    y=yo2+fd*B*sin(course-1.5*pi);
%     x=xo2+fd*B*sin(course-0.5*pi);
%     y=yo2+fd*B*cos(course-0.5*pi);
    plt5 = plot(x,y,'-b');
    hold on
    xo3=xo-0.5*B*sin(course)-0.5*L*cos(course);
    yo3=yo-0.5*L*sin(course)+0.5*B*cos(course);
    x=xo3+fd*0.75*L*cos(course);
    y=yo3+fd*0.75*L*sin(course);
    plt6 = plot(x,y,'-b');
    axis equal
    hold on
    







