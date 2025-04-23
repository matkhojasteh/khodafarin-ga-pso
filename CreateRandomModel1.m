function model1=CreateRandomModel1()

load('var.mat','b_canal','h_canal','H_innerslope','V_innerslope',...
    'W_leftberm','W_rightberm','H_total','Hmin_fillingsteps','Wmin_berms');

   
    XC=499.25;
    XLC=[497,489.5,484.5];
    YLC=[25,30,30];
    XRC=[503,510.5,514];
    YRC=YLC;
    
    phi_max=atan(1);
    phi_min=atan(1/3);

% section
    
    for l=1:3
        XL(l)=XLC(l);
        YL(l)=YLC(l);
        XR(l)=XRC(l);
        YR(l)=YRC(l);
    end
    
    %% for h1>10
    h1=H_total;
    
        I=numel(XLC);
        e=h1/tan((phi_max-phi_min)*rand() + phi_min);
        XL(I+1)=XL(I)-e;
        YL(I+1)=YL(I)-h1;
        
        zX=XL(4);
        zY=YL(4);   
    
  N=numel(zX);
   
   model1.N=N;            
   model1.Position.zX=zX;
   model1.Position.zY=zY;
   model1.xmin=XL(3)-h1/tan(phi_min);
   model1.xmax=XL(3)-h1/tan(phi_max);
   model1.ymin=h1;
   model1.ymax=h1;
   
end