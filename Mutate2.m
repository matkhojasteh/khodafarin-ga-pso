function y=Mutate2(x,mu)


   load('var.mat','Wmin_berms');
        
    nVar=numel(x.Position.zX);
    
    nmu=ceil(mu*nVar);
    
    j=randsample(nVar,nmu);
    
     x=[x.Position.zX;x.Position.zY];
    
    max1=484.5-(5)/tan(1.107);
    max2=max1-2;
    max3=max2-(5)/tan(1.107);
    max4=max3-2;
    max5=max4-(5)/tan(1.107);
        VarMin.x=[max1 max2 max3 max4 max5];     
    min1=484.5-(5)/tan(.1974);
    min2=min1-2;
    min3=min2-(5)/tan(.1974); 
    min4=max3-2;
    min5=max4-(5)/tan(.1974);
        VarMax.x=[min1 min2 min3 min4 min5];           
     
          
     sigma_X=0.1.*(VarMax.x-VarMin.x);
          
     
    t=x;
    t(1,j)=x(1,j)+sigma_X(1,j)*randn(size(j));
    
    
    t(1,:)=max(t(1,:),VarMin.x);
    t(1,:)=min(t(1,:),VarMax.x);
   
    t(1,:)=sort(t(1,:),'descend');
   if t(1,1)-t(1,2)<2
        t(1,2)=t(1,1)-2;
   end
    if t(1,3)-t(1,4)<2
        t(1,4)=t(1,3)-2;
    end
    if t(1,2)-t(1,3)<4.52
        t(1,3)=t(1,2)-4.52;
    end
    if t(1,4)-t(1,5)<4.52
        t(1,5)=t(1,4)-4.52;
    end
    
    y=[ ];
    
    
    y.zX=t(1,:);
    y.zY=t(2,:);
    
   
end
    
     
      
  
    
    