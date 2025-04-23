function [Cost sol]= MyCost(f)

   global NFE;
   if isempty(NFE)
      NFE=0;
   end

   NFE=NFE+1;
    
    XC=499.25;
    XLC=[497,489.5,484.5];
    YLC=[25,30,30];
    XRC=[503,510.5,514];
    YRC=YLC;
    
    zX=f.Position.zX;
    zY=f.Position.zY;
    
    e=numel(zX);
    
    XL=[XLC zX];
    YL=[YLC zY];
    
    XR=XRC;
        
    for i=1:e
    
        
    XR(i+3)=XL(i+3)+2*(XC-XL(i+3));
    
    end
    
    YR=YL;
    
    x=[XL fliplr(XR) XL(1)];
    y=[YL fliplr(YR) YL(1)];
    
   
    n = numel(x);

   
    Cost = 0;
    for k = 1:n-1
               
        Cost = Cost + x(k)*y(k+1) - y(k)*x(k+1);
        
    end
    Cost = Cost/2;
    sol.x=x;    
    sol.y=y;
    
end