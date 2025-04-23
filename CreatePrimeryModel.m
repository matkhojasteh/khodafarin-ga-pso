function pmodel=CreatePrimeryModel()

load('var.mat','b_canal','h_canal','H_innerslope','V_innerslope',...
    'W_leftberm','W_rightberm','H_total','Hmin_fillingsteps','Wmin_berms');

%Constants

Alpha=atan(V_innerslope/H_innerslope);

% center
XC_canal=500;
XC=499.25;
YC=H_total+(H_total-h_canal);

%left side
XL(1)=XC_canal-b_canal/2;
YL(1)=YC;

XL(2)=XL(1)-h_canal/tan(Alpha);
YL(2)=YL(1)+h_canal;

XL(3)=XL(2)-W_leftberm;
YL(3)=YL(2);

%right side
XR(1)=XC_canal+b_canal/2;
YR(1)=YC;

XR(2)=XR(1)+h_canal/tan(Alpha);
YR(2)=YR(1)+h_canal;

XR(3)=XR(2)+W_rightberm;
YR(3)=YR(2);

%primery section
for I=4:8
   if mod(I,2)==0
       XL(I)=XL(I-1)-h_canal/tan(Alpha);
       YL(I)=YL(I-1)-h_canal;
       
       XR(I)=XL(I)+2*(XC-XL(I));
       YR(I)=YL(I);
   else
       XL(I)=XL(I-1)-Wmin_berms;
       YL(I)=YL(I-1);
       
       XR(I)=XL(I)+2*(XC-XL(I));
       YR(I)=YL(I);   
   end
end

xsoil1=[XL fliplr(XR) XL(1)];
ysoil1=[YL fliplr(YR) YL(1)];

XLC=XL(1:3);
YLC=YL(1:3);

XRC=XR(1:3);
YRC=YR(1:3);

zX=XL(4:end);
zY=YL(4:end);


    N=numel(zX);
   
    pmodel.N=N;
    pmodel.Position.zX=zX;
    pmodel.Position.zY=zY;
    
end