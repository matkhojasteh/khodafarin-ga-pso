function model2=CreateRandomModel2()

load('var.mat','b_canal','h_canal','H_innerslope','V_innerslope',...
    'W_leftberm','W_rightberm','H_total','Hmin_fillingsteps',...
    'Wmin_berms');

    XC=499.25;
    XLC=[497,489.5,484.5];
    YLC=[25,30,30];
    XRC=[503,510.5,514];
    YRC=YLC;

    phi_max=atan(1/1.4);
    phi_min=atan(1/3);

% section
    
    for l=1:3
        XL(l)=XLC(l);
        YL(l)=YLC(l);
        XR(l)=XRC(l);
        YR(l)=YRC(l);
    end
    
    %% for h1<10
    h1=randi([6 9],1);
    h2=H_total-h1 ;
        I=numel(XLC);
        e1=h1/tan((phi_max-phi_min)*rand() + phi_min);
        XL(I+1)=XL(I)-e1;
        YL(I+1)=YL(I)-h1;
        w = round((10-Wmin_berms)*rand()+Wmin_berms);
        XL(I+2)=XL(I+1)-w;
        YL(I+2)=YL(I+1);
        e2=h2/tan((phi_max-phi_min)*rand() + phi_min);
        XL(I+3)=XL(I+2)-e2;
        YL(I+3)=YL(I+2)-h2;
        
        zX=XL(4:end);
        zY=YL(4:end);   
        Vmax=6/tan( phi_max);
        Vmin=9/tan( phi_min);
               
  N=numel(zX);
   
   model2.N=N;               
   model2.Position.zX=zX;
   model2.Position.zY=zY;     
   model2.x=XL(3);  
   model2.Vmax=Vmax;
   model2.Vmin=Vmin;
        
end