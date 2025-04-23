function [y1 y2]=SinglePointCrossover2(x1,x2)
    
     
    x1=[x1.Position.zX;x1.Position.zY];
    x2=[x2.Position.zX;x2.Position.zY];
    
    x11=x1(1,:);
    x12=x1(2,:);
    x21=x2(1,:);
    x22=x2(2,:);
    
    if x11(1,3)-2<x21(1,4)
    x21(1,4)=x11(1,3)-2;
    x21(1,5)=x21(1,4)-(x22(1,4)-x22(1,5))/tan(1.107);
   end
   
    t1=[x11(1:3) x21(4:5);x12(1:5)];
    
    t2=[x21(1:3) x11(4:5);x22(1:5)];

    y1=[];
    y2=[];
    
    y1.zX=t1(1,:);
    y1.zY=t1(2,:);
       
    y2.zX=t2(1,:);
    y2.zY=t2(2,:);
    
end