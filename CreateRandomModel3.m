function model3=CreateRandomModel3()

load('var.mat','b_canal','h_canal','H_innerslope','V_innerslope',...
    'W_leftberm','W_rightberm','H_total','Hmin_fillingsteps',...
    'Wmin_berms');

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
    
    %% for h1=5 && h2=5
    h=5;
    
    
        I=numel(XLC);
        e1=h/tan((phi_max-phi_min)*rand() + phi_min);
        XL(I+1)=XL(I)-e1;
        YL(I+1)=YL(I)-h;
        w1 = round((10-Wmin_berms)*rand()+Wmin_berms);
        XL(I+2)=XL(I+1)-w1;
        YL(I+2)=YL(I+1);
        e2=h/tan((phi_max-phi_min)*rand() + phi_min);
        XL(I+3)=XL(I+2)-e2;
        YL(I+3)=YL(I+2)-h;
         w2 = round((10-Wmin_berms)*rand()+Wmin_berms);
        XL(I+4)=XL(I+3)-w2;
        YL(I+4)=YL(I+3);
        e3=h/tan((phi_max-phi_min)*rand() + phi_min);
        XL(I+5)=XL(I+4)-e3;
        YL(I+5)=YL(I+4)-h;        
       
        
        zX=XL(4:end);
        zY=YL(4:end);   
        Vmax=5/tan( phi_max);
        Vmin=5/tan( phi_min);
               
    N=numel(zX);
   
   model3.N=N;               
   model3.Position.zX=zX;
   model3.Position.zY=zY;     
   model3.x=XL(3);
    model3.Vmax=Vmax;
   model3.Vmin=Vmin;
        
end