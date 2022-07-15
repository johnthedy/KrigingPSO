function g=G(x,probtype)
if probtype==1
    k=4;
    m=7;
    g=min(min(0.1*(x(1)-x(2))^2-(x(1)+x(2))/sqrt(2)+k,0.1*(x(1)-x(2))^2+(x(1)+x(2))/sqrt(2)+k),min(x(1)-x(2)+m/sqrt(2),-x(1)+x(2)+m/sqrt(2)));
elseif probtype==2
    k=5;
    m=8;
    g=min(min(0.1*(x(1)-x(2))^2-(x(1)+x(2))/sqrt(2)+k,0.1*(x(1)-x(2))^2+(x(1)+x(2))/sqrt(2)+k),min(x(1)-x(2)+m/sqrt(2),-x(1)+x(2)+m/sqrt(2)));
elseif probtype==3
    k=5;
    m=9;
    g=min(min(0.1*(x(1)-x(2))^2-(x(1)+x(2))/sqrt(2)+k,0.1*(x(1)-x(2))^2+(x(1)+x(2))/sqrt(2)+k),min(x(1)-x(2)+m/sqrt(2),-x(1)+x(2)+m/sqrt(2)));
elseif probtype==4
    g=min(min(0.1*(x(1)-x(2))^2-(x(1)+x(2))/sqrt(2)+3,0.1*(x(1)-x(2))^2+(x(1)+x(2))/sqrt(2)+3),min(x(1)-x(2)+3.5*sqrt(2),-x(1)+x(2)+3.5*sqrt(2)));
elseif probtype==5
    g=0.5*(x(1)-2)^2-1.5*(x(2)-5)^3-3;
elseif probtype==6
    g=2.5-0.2357*(x(1)-x(2))+0.00463*(x(1)+x(2)-20)^4;
elseif probtype==7
    g=-0.5*(x(1)-x(2))^2-(x(1)+x(2))/sqrt(2)+3;
elseif probtype==8
    w=sqrt((x(2)+x(3))/x(1));
    g=3*x(4)-abs(2*x(6)/(x(1)*w^2)*sin(0.5*w*x(5)));
elseif probtype==9
    w=sqrt((x(2)+x(3))/x(1));
    g=3*x(4)-abs(2*x(6)/(x(1)*w^2)*sin(0.5*w*x(5)));
elseif probtype==10
    w=sqrt((x(2)+x(3))/x(1));
    g=3*x(4)-abs(2*x(6)/(x(1)*w^2)*sin(0.5*w*x(5)));
elseif probtype==11
    dummy1=2*x(1)+2*x(3)-4.5*x(4);
    dummy2=2*x(1)+x(2)+x(3)-4.5*x(4);
    dummy3=x(1)+x(2)+2*x(3)-4.5*x(4);
    dummy4=x(1)+2*x(2)+x(3)-4.5*x(4);
    g=min([dummy1,dummy2,dummy3,dummy4]);
elseif probtype==12
    g=0.03-(x(1)*(x(2)^2)/2)*(3.81/(x(3)*x(4))+1.13/(x(5)*x(6)));
elseif probtype==13
    dummy0=0;
    for i=1:9
        dummy0=dummy0+(x(i)^2);
    end
    g=2+0.015*dummy0-x(10);
elseif probtype==14
    x(10)=deg2rad(5);x(11)=deg2rad(5);
    x(5)=x(5)*1000;
    x(6)=x(6)*1000;
    x(7)=x(7)*1000;
    x(8)=x(8)*1000;
    
    M=x(5)*x(3)*cos(x(10))+x(6)*x(4)*cos(x(11));
    A=(pi/4)*(x(2)^2-(x(2)-2*x(1))^2);
    I=(pi/64)*(x(2)^4-(x(2)-2*x(1))^4);J=2*I;
    
    szx=(x(8)*x(2))/(2*J);
    sx=(x(7)+x(5)*sin(x(10))+x(6)*sin(x(11)))/A+(M*x(2))/(2*I);
    g=x(9)-sqrt((sx)^2+3*(szx)^2);
end
end