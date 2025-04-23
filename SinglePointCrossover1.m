function [y1 y2]=SinglePointCrossover1(x1,x2)
    
     
    x1=[x1.Position.zX;x1.Position.zY];
    x2=[x2.Position.zX;x2.Position.zY];
    
    x11=x1(1,:);
    x12=x1(2,:);
    x21=x2(1,:);
    x22=x2(2,:);
    
   if x11(1,1)-2<x21(1,2)
    x21(1,2)=x11(1,1)-2;
    x21(1,3)=x21(1,2)-(x22(1,2)-x22(1,3))/tan(1.107);
   end
   t1=[x11(1:1) x21(2:3);x12(1:3)];
    
    t2=[x21(1:1) x11(2:3);x22(1:3)];

    y1=[];
    y2=[];
    
    y1.zX=t1(1,:);
    y1.zY=t1(2,:);
       
    y2.zX=t2(1,:);
    y2.zY=t2(2,:);
    
end