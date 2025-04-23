clc;
clear;
close all;

% Input Datas:
prompt = {'Enter b_canal(m):','Enter h_canal(m):','Enter H_innerslope:','Enter V_innerslope:',...
    'Enter W_leftberm(m):','Enter W_rightberm(m):','Enter H_total(m):',...
    'Enter Hmin_fillingsteps(m):','Enter Wmin_berms(m):'};
    
title = 'Input Datas';
lines = 1;
def={'6','5','1.5','1','5','3.5','15','5','2'};
answer = inputdlg(prompt, title, lines,def);
assignin('base', 'b_canal', answer{1});
b_canal=str2num(b_canal);
assignin('base', 'h_canal', answer{2});
h_canal=str2num(h_canal);
assignin('base', 'H_innerslope', answer{3});
H_innerslope=str2num(H_innerslope);
assignin('base', 'V_innerslope', answer{4});
V_innerslope=str2num(V_innerslope);
assignin('base', 'W_leftberm', answer{5});
W_leftberm=str2num(W_leftberm);
assignin('base', 'W_rightberm', answer{6});
W_rightberm=str2num(W_rightberm);
assignin('base', 'H_total', answer{7});
H_total=str2num(H_total);
assignin('base', 'Hmin_fillingsteps', answer{8});
Hmin_fillingsteps=str2num(Hmin_fillingsteps);
assignin('base', 'Wmin_berms', answer{9});
Wmin_berms=str2num(Wmin_berms);
  


save('var.mat','b_canal','h_canal','H_innerslope','V_innerslope',...
    'W_leftberm','W_rightberm','H_total','Hmin_fillingsteps','Wmin_berms');


%% Problem Definition

global NFE;
NFE=0;

pmodel=CreatePrimeryModel();
model1=CreateRandomModel1();
model2=CreateRandomModel2();
model3=CreateRandomModel3();

CostFunction=@(f) MyCost(f);    % Cost Function

%% PSO_GA Parameters

MaxIt=2;      % Maximum Number of Iterations

nPop=10000;        % Population Size (Swarm Size)

phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;            % Inertia Weight
wdamp=0.999;     % Inertia Weight Damping Ratio
c1=phi1*chi;           % Personal Learning Coefficient
c2=phi2*chi;           % Global Learning Coefficient

pc=0.3;                 % Crossover Percentage
nc=2*round(pc*nPop/2);  % Number of Offsprings (Parnets)

pm=0.1;                 % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants

mu=0.05;         % Mutation Rate

ANSWER=questdlg('Choose selection method:','Genetic Algorith',...
    'Roulette Wheel','Tournament','Random','Roulette Wheel');

UseRouletteWheelSelection=strcmp(ANSWER,'Roulette Wheel');
UseTournamentSelection=strcmp(ANSWER,'Tournament');
UseRandomSelection=strcmp(ANSWER,'Random');

if UseRouletteWheelSelection
    beta=10;         % Selection Pressure
end

if UseTournamentSelection
    TournamentSize=3;   % Tournamnet Size
end

pause(0.1);

%% Initialization

empty_pop.Position=[ ];
empty_pop.Cost=[ ];
empty_pop.sol=[ ];
empty_pop.Velocity=[ ];
empty_pop.Best.Position=[ ];
empty_pop.Best.Cost=[ ];
empty_pop.Best.sol=[ ];

pop1=repmat(empty_pop,nPop,1);
pop2=repmat(empty_pop,nPop,1);
pop3=repmat(empty_pop,nPop,1);

[pmodel.cost pmodel.sol]=MyCost(pmodel);

GlobalBest.Cost=pmodel.cost;
GlobalBest1.Cost=inf;
GlobalBest2.Cost=inf;
GlobalBest3.Cost=inf;

%% for h1>10

nVar=model1.N;       % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix
VarMin.x=model1.xmin;           % Lower Bound of Variables
VarMax.x=model1.xmax;           % Upper Bound of Variables
VarMin.y=model1.ymin;           % Lower Bound of Variables
VarMax.y=model1.ymax;           % Upper Bound of Variables

% Velocity Limits
alpha=0.01;
VelMax.x=alpha*(VarMax.x-VarMin.x);    % Maximum Velocity
VelMin.x=-VelMax.x;                    % Minimum Velocity
VelMax.y=alpha*(VarMax.y-VarMin.y);    % Maximum Velocity
VelMin.y=-VelMax.y;                    % Minimum Velocity

for i=1:500
    
    % Initialize Position
    pop1(i).Position.zX=unifrnd(VarMin.x,VarMax.x,VarSize);
    pop1(i).Position.zY=unifrnd(VarMin.y,VarMax.y,VarSize);
    
    % Initialize Velocity
    pop1(i).Velocity.zX=zeros(VarSize);
    pop1(i).Velocity.zY=zeros(VarSize);
    
    % Evaluation
    [pop1(i).Cost pop1(i).sol]=MyCost(pop1(i));
    
    % Update Personal Best
    pop1(i).Best.Position=pop1(i).Position;
    pop1(i).Best.Cost=pop1(i).Cost;
    pop1(i).Best.sol=pop1(i).sol;
    
    % Update Global Best
    if pop1(i).Best.Cost<GlobalBest1.Cost
        
        GlobalBest1=pop1(i).Best;
        
    end
    
end

BestCost1=zeros(MaxIt,1);

nfe=zeros(MaxIt,1);


%% PSO_GA Main Loop

for it=1:10
    %PSO Operetors
    for i=1:500
        
      % x Part
        
        % Update Velocity
        pop1(i).Velocity.zX = w*pop1(i).Velocity.zX ...
            + c1*rand(VarSize).*(pop1(i).Best.Position.zX-pop1(i).Position.zX) ...
            + c2*rand(VarSize).*(GlobalBest1.Position.zX-pop1(i).Position.zX);
        
        % Update Velocity Bounds
        pop1(i).Velocity.zX = max(pop1(i).Velocity.zX,VelMin.x);
        pop1(i).Velocity.zX = min(pop1(i).Velocity.zX,VelMax.x);
        
        % Update Position
        pop1(i).Position.zX = pop1(i).Position.zX + pop1(i).Velocity.zX;
        
        % Velocity Mirroring
        OutOfTheRange=(pop1(i).Position.zX<VarMin.x | pop1(i).Position.zX>VarMax.x);
        pop1(i).Velocity.zX(OutOfTheRange)=-pop1(i).Velocity.zX(OutOfTheRange);
        
        % Update Position Bounds
        pop1(i).Position.zX = max(pop1(i).Position.zX,VarMin.x);
        pop1(i).Position.zX = min(pop1(i).Position.zX,VarMax.x);
        
        % y Part
        
        % Update Velocity
        pop1(i).Velocity.zY = w*pop1(i).Velocity.zY ...
            + c1*rand(VarSize).*(pop1(i).Best.Position.zY-pop1(i).Position.zY) ...
            + c2*rand(VarSize).*(GlobalBest1.Position.zY-pop1(i).Position.zY);
        
        % Update Velocity Bounds
        pop1(i).Velocity.y = max(pop1(i).Velocity.zY,VelMin.y);
        pop1(i).Velocity.y = min(pop1(i).Velocity.zY,VelMax.y);
        
        % Update Position
        pop1(i).Position.zY = pop1(i).Position.zY + pop1(i).Velocity.zY;
        
        % Velocity Mirroring
        OutOfTheRange=(pop1(i).Position.zY<VarMin.y | pop1(i).Position.zY>VarMax.y);
        pop1(i).Velocity.zY(OutOfTheRange)=-pop1(i).Velocity.zY(OutOfTheRange);
        
        % Update Position Bounds
        pop1(i).Position.zY = max(pop1(i).Position.zY,VarMin.y);
        pop1(i).Position.zY = min(pop1(i).Position.zY,VarMax.y);
        
        % Evaluation
       [pop1(i).Cost pop1(i).sol] =MyCost(pop1(i));
        
        % Update Personal Best
        if pop1(i).Cost<pop1(i).Best.Cost
            
            pop1(i).Best.Position=pop1(i).Position;
            pop1(i).Best.Cost=pop1(i).Cost;
            pop1(i).Best.sol=pop1(i).sol;    
            
            % Update Global Best
            if pop1(i).Best.Cost<GlobalBest1.Cost
                GlobalBest1=pop1(i).Best;
            end
            
        end
        
        
    end
     [~, order]=sort([pop1(:).Cost]);
   pop1=pop1(order);
    BestCost1(it)=GlobalBest1.Cost;  
          
    nfe(it)=NFE;
    
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost1 = ' num2str(BestCost1(it))]);
    
    w=w*wdamp;
    
 end

%% Results

figure;
plot(nfe,BestCost1,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost1');

%% for h1<10

nVar=model2.N;       % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix
Var.x=model2.x;         
Max=model2.Vmax;
Min=model2.Vmin;
a1=Var.x-Max;
a2=a1-2;
a3=a2-Max;
b1=Var.x-Min;
b2=b1-10;
b3=b2-Min;
VarMax.x=[a1 a2 a3];
VarMin.x=[b1 b2 b3];

% Velocity Limits
alpha=0.1;
r1=alpha*(a1-b1);
r2=alpha*((r1-2)-(r1-10));
r3=alpha*((r2-Max)-(r2-Min));
VelMax.x=[r1 r2 r3];                     % Maximum Velocity
VelMin.x=-VelMax.x;                    % Minimum Velocity
 q=Min;
qq=10;
for i=1:nPop
    
    % Initialize Position
    o1=Var.x-q;
    o2=o1-qq;
    o3=o2-q;
    pop2(i).Position.zX=[o1 o2 o3];
    a=randi([21 24],1);
    pop2(i).Position.zY=[a a 15];
    % Initialize Velocity
    pop2(i).Velocity.zX=zeros(VarSize);
    pop2(i).Velocity.zY=zeros(VarSize);
    
    % Evaluation
    pop2(i).Cost=MyCost(pop2(i));
    
    % Update Personal Best
    pop2(i).Best.Position=pop2(i).Position;
    pop2(i).Best.Cost=pop2(i).Cost;
       
    % Update Global Best
    if pop2(i).Best.Cost<GlobalBest2.Cost
        
        GlobalBest2=pop2(i).Best;
        
    end
     q=q-(Min-Max)/nPop;
    qq=qq-(10-2)/nPop;
end

BestCost2=zeros(MaxIt,1);

nfe=zeros(MaxIt,1);


%% PSO_GA Main Loop

for it=1:MaxIt
    
    %PSO Operetors
    for i=1:nPop
        
      % x Part
        
        % Update Velocity
        pop2(i).Velocity.zX = w*pop2(i).Velocity.zX ...
            + c1*rand(VarSize).*(pop2(i).Best.Position.zX-pop2(i).Position.zX) ...
            + c2*rand(VarSize).*(GlobalBest2.Position.zX-pop2(i).Position.zX);
        
        % Update Velocity Bounds
        pop2(i).Velocity.zX = max(pop2(i).Velocity.zX,VelMin.x);
        pop2(i).Velocity.zX = min(pop2(i).Velocity.zX,VelMax.x);
        
        % Update Position
        pop2(i).Position.zX = pop2(i).Position.zX + pop2(i).Velocity.zX;
        
        % Velocity Mirroring
        OutOfTheRange=(pop2(i).Position.zX<VarMin.x | pop2(i).Position.zX>VarMax.x);
        pop2(i).Velocity.zX(OutOfTheRange)=-pop2(i).Velocity.zX(OutOfTheRange);
        
        % Update Position Bounds
        pop2(i).Position.zX = max(pop2(i).Position.zX,VarMin.x);
        pop2(i).Position.zX = min(pop2(i).Position.zX,VarMax.x);
        pop2(i).Position.zX = sort(pop2(i).Position.zX,'descend');
%          if pop2(i).Position.zX (1,2)-pop2(i).Position.zX (1,1)<2 
%            pop2(i).Position.zX (1,2)= pop2(i).Position.zX (1,1)-2;
%         end
%         if pop2(i).Position.zX (1,3)-pop2(i).Position.zX (1,2)<Max
%            pop2(i).Position.zX (1,3)= pop2(i).Position.zX (1,2)-Max;
%         end

        
        % Evaluation
       pop2(i).Cost =MyCost(pop2(i));
        
        % Update Personal Best
        if pop2(i).Cost<pop2(i).Best.Cost
            
            pop2(i).Best.Position=pop2(i).Position;
            pop2(i).Best.Cost=pop2(i).Cost;
                       
            % Update Global Best
            if pop2(i).Best.Cost<GlobalBest2.Cost
                GlobalBest2=pop2(i).Best;
            end
            
        end
        
        
    end
     
    %GA Operators
    % Sort Population
Costs2=[pop2.Cost];
[Costs2, SortOrder]=sort(Costs2);
pop2=pop2(SortOrder);

% Store Best Solution
BestSol2=pop2(1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Store Cost
WorstCost=pop2(end).Cost;

      % Calculate Selection Probabilities
    if UseRouletteWheelSelection
        P=exp(-beta*Costs2/WorstCost);
        P=P/sum(P);
    end
    
    % Crossover
    popc=repmat(empty_pop,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices
        if UseRouletteWheelSelection
            i1=RouletteWheelSelection(P);
            i2=RouletteWheelSelection(P);
        end
        if UseTournamentSelection
            i1=TournamentSelection(pop2,TournamentSize);
            i2=TournamentSelection(pop2,TournamentSize);
        end
        if UseRandomSelection
            i1=randi([1 nPop]);
            i2=randi([1 nPop]);
        end

        % Select Parents
        p1=pop2(i1);
        p2=pop2(i2);
        
        % Apply Crossover
        [popc(k,1).Position popc(k,2).Position]=SinglePointCrossover1(p1,p2);
               
        % Evaluate Offsprings
        [popc(k,1).Cost popc(k,1).sol] =MyCost(popc(k,1));
       
        [popc(k,2).Cost popc(k,2).sol]=MyCost(popc(k,2));
        
        if p1.Best.Cost<p2.Best.Cost
            popc(k,1).Best=p1.Best;
            popc(k,2).Best=p1.Best;
        else
            popc(k,1).Best=p2.Best;
            popc(k,2).Best=p2.Best;
        end
        
        if rand<.5
            popc(k,1).Velocity=p1.Velocity;
            popc(k,2).Velocity=p2.Velocity;
        else
            popc(k,1).Velocity=p2.Velocity;
            popc(k,2).Velocity=p1.Velocity;
      
        end
    end
    popc=popc(:);
    
    
    % Mutation
    popm=repmat(empty_pop,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop2(i);
        
        % Apply Mutation
        popm(k).Position=Mutate1(p,mu);
        
        % Evaluate Mutant
        [popm(k).Cost popm(k).sol]=MyCost(popm(k));
        popm(k).Velocity=p.Velocity;
        popm(k).Best=p.Best;
        
    end
    
    % Create Merged Population
    pop2=[pop2
         popc
         popm];
     
    % Sort Population
    Costs2=[pop2.Cost];
    [Costs2, SortOrder]=sort(Costs2);
    pop2=pop2(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop2(end).Cost);
    
    % Truncation
    pop2=pop2(1:nPop);
  
       % Update Personal Best
        if pop2(i).Cost<pop2(i).Best.Cost
            
            pop2(i).Best.Position=pop2(i).Position;
            pop2(i).Best.Cost=pop2(i).Cost;
                       
            % Update Global Best
            if pop2(i).Best.Cost<GlobalBest2.Cost
                GlobalBest2=pop2(i).Best;
            end
            
        end 
     [~, order]=sort([pop2(:).Cost]);
    pop2=pop2(order);   
    BestCost2(it)=GlobalBest2.Cost;
    
    nfe(it)=NFE;
    
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost2 = ' num2str(BestCost2(it))]);
    
    w=w*wdamp;
    
       
end

%% Results

figure;
plot(nfe,BestCost2,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost2');

%% for h1=5 && h2=5

nVar=model3.N;       % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix
Var.x=model3.x;         
Max=model3.Vmax;
Min=model3.Vmin;
a1=Var.x-Max;
a2=a1-2;
a3=a2-Max;
a4=a3-2;
a5=a4-Max;
b1=Var.x-Min;
b2=b1-10;
b3=b2-Min;
b4=b3-10;
b5=b4-Min;
VarMax.x=[a1 a2 a3 a4 a5];
VarMin.x=[b1 b2 b3 b4 b5];


% Velocity Limits
alpha=0.1;
r1=alpha*(a1-b1);
r2=alpha*((r1-2)-(r1-10));
r3=alpha*((r2-Max)-(r2-Min));
r4=alpha*((r3-2)-(r3-10));
r5=alpha*((r4-Max)-(r4-Min));
VelMax.x=[r1 r2 r3 r4 r5];                        % Maximum Velocity
VelMin.x=-VelMax.x;                       % Minimum Velocity
q=Min;
qq=10;
for i=1:nPop
    
    % Initialize Position
    o1=Var.x-q;
    o2=o1-qq;
    o3=o2-q;
    o4=o3-qq;
    o5=o4-q;
    pop3(i).Position.zX=[o1 o2 o3 o4 o5];
    pop3(i).Position.zY=[25 25 20 20 15];
    
    % Initialize Velocity
    pop3(i).Velocity.zX=zeros(VarSize);
    pop3(i).Velocity.zY=zeros(VarSize);
    
    % Evaluation
    pop3(i).Cost=MyCost(pop3(i));
    
    % Update Personal Best
    pop3(i).Best.Position=pop3(i).Position;
    pop3(i).Best.Cost=pop3(i).Cost;
       
    % Update Global Best
    if pop3(i).Best.Cost<GlobalBest3.Cost
        
        GlobalBest3=pop3(i).Best;
        
    end
    q=q-(Min-Max)/nPop;
    qq=qq-(10-2)/nPop;
end

BestCost3=zeros(MaxIt,1);

nfe=zeros(MaxIt,1);


%% PSO_GA Main Loop

for it=1:MaxIt
    %PSO Operetors
    for i=1:nPop
        
      % x Part
        
        % Update Velocity
        pop3(i).Velocity.zX = w*pop3(i).Velocity.zX ...
            + c1*sort(rand(VarSize)).*(pop3(i).Best.Position.zX-pop3(i).Position.zX) ...
            + c2*sort(rand(VarSize)).*(GlobalBest3.Position.zX-pop3(i).Position.zX);
        
        % Update Velocity Bounds
        pop3(i).Velocity.zX = max(pop3(i).Velocity.zX,VelMin.x);
        pop3(i).Velocity.zX = min(pop3(i).Velocity.zX,VelMax.x);
        
        % Update Position
        pop3(i).Position.zX = pop3(i).Position.zX + pop3(i).Velocity.zX;
        
        % Velocity Mirroring
        OutOfTheRange=(pop3(i).Position.zX<VarMin.x | pop3(i).Position.zX>VarMax.x);
        pop3(i).Velocity.zX(OutOfTheRange)=-pop3(i).Velocity.zX(OutOfTheRange);
        
        % Update Position Bounds
        pop3(i).Position.zX = max(pop3(i).Position.zX,VarMin.x);
        pop3(i).Position.zX = min(pop3(i).Position.zX,VarMax.x);
        pop3(i).Position.zX = sort(pop3(i).Position.zX,'descend');
%          if pop3(i).Position.zX (1,2)-pop3(i).Position.zX (1,1)<2 
%            pop3(i).Position.zX (1,2)= pop3(i).Position.zX (1,1)-2;
%         end
%         if pop3(i).Position.zX (1,3)-pop3(i).Position.zX (1,2)<Max
%            pop3(i).Position.zX (1,3)= pop3(i).Position.zX (1,2)-Max;
%         end
%         if pop3(i).Position.zX (1,4)-pop3(i).Position.zX (1,3)<2
%            pop3(i).Position.zX (1,4)= pop3(i).Position.zX (1,3)-2;
%         end
%         if pop3(i).Position.zX (1,5)-pop3(i).Position.zX (1,4)<Max
%            pop3(i).Position.zX (1,5)= pop3(i).Position.zX (1,4)-Max;
%         end

        
        % Evaluation
       pop3(i).Cost =MyCost(pop3(i));
        
        % Update Personal Best
        if pop3(i).Cost<pop3(i).Best.Cost
            
            pop3(i).Best.Position=pop3(i).Position;
            pop3(i).Best.Cost=pop3(i).Cost;
                       
            % Update Global Best
            if pop3(i).Best.Cost<GlobalBest3.Cost
                GlobalBest3=pop3(i).Best;
            end
            
        end
        
        
    end
    
    %GA Operators
      % Sort Population
    Costs3=[pop3.Cost];
    [Costs3, SortOrder]=sort(Costs3);
    pop3=pop3(SortOrder);

    % Store Best Solution
    BestSol3=pop3(1);

    % Array to Hold Best Cost Values
    BestCost=zeros(MaxIt,1);

    % Store Cost
    WorstCost=pop3(end).Cost;

         % Calculate Selection Probabilities
    if UseRouletteWheelSelection
        P=exp(-beta*Costs3/WorstCost);
        P=P/sum(P);
    end
    
    % Crossover
    popc=repmat(empty_pop,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices
        if UseRouletteWheelSelection
            i1=RouletteWheelSelection(P);
            i2=RouletteWheelSelection(P);
        end
        if UseTournamentSelection
            i1=TournamentSelection(pop3,TournamentSize);
            i2=TournamentSelection(pop3,TournamentSize);
        end
        if UseRandomSelection
            i1=randi([1 nPop]);
            i2=randi([1 nPop]);
        end

        % Select Parents
        p1=pop3(i1);
        p2=pop3(i2);
        
        % Apply Crossover
        [popc(k,1).Position popc(k,2).Position]=SinglePointCrossover2(p1,p2);
               
        % Evaluate Offsprings
        [popc(k,1).Cost popc(k,1).sol]=MyCost(popc(k,1));
        
        [popc(k,2).Cost popc(k,2).sol]=MyCost(popc(k,2));
       
          if p1.Best.Cost<p2.Best.Cost
            popc(k,1).Best=p1.Best;
            popc(k,2).Best=p1.Best;
         else
            popc(k,1).Best=p2.Best;
            popc(k,2).Best=p2.Best;
         end
        
        if rand<.5
            popc(k,1).Velocity=p1.Velocity;
            popc(k,2).Velocity=p2.Velocity;
        else
            popc(k,1).Velocity=p2.Velocity;
            popc(k,2).Velocity=p1.Velocity;
        end
 
    end
    popc=popc(:);
    
    
    % Mutation
    popm=repmat(empty_pop,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop3(i);
        
        % Apply Mutation
        popm(k).Position=Mutate2(p,mu);
        
        % Evaluate Mutant
        [popm(k).Cost popm(k).sol]=MyCost(popm(k));
        popm(k).Velocity=p.Velocity;
        popm(k).Best=p.Best;
    end
    
    % Create Merged Population
    pop3=[pop3
         popc
         popm];
     
    % Sort Population
    Costs3=[pop3.Cost];
    [Costs3, SortOrder]=sort(Costs3);
    pop3=pop3(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop3(end).Cost);
    
    % Truncation
    pop3=pop3(1:nPop);
      
        % Update Personal Best
        if pop3(i).Cost<pop3(i).Best.Cost
            
            pop3(i).Best.Position=pop3(i).Position;
            pop3(i).Best.Cost=pop3(i).Cost;
                       
            % Update Global Best
            if pop3(i).Best.Cost<GlobalBest3.Cost
                GlobalBest3=pop3(i).Best;
            end
            
        end 
    
     [~, order]=sort([pop3(:).Cost]);
    pop3=pop3(order);   
    BestCost3(it)=GlobalBest3.Cost;
    
    nfe(it)=NFE;
    
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it)) ', Best Cost3 = ' num2str(BestCost3(it))]);
    
    w=w*wdamp;
    
      
end

%% Results

figure;
plot(nfe,BestCost3,'LineWidth',2);
xlabel('NFE');
ylabel('Best Cost3');

n=numel(pop1);
P1=struct2cell(pop1);
for i=1:n
   PP1{i}=struct2cell(P1{1,i});
   Pop1{i}=[PP1{1,i}{1,1} 
                  PP1{1,i}{2,1}];
    C1{i}=P1{2,i};           
end
Pop1=(cell2mat(Pop1))';
Cost1=(cell2mat(C1))';

n=numel(pop2);
P2=struct2cell(pop2);
for i=1:n
   PP2{i}=struct2cell(P2{1,i});
   Pop2{i}=[PP2{1,i}{1,1} 
                  PP2{1,i}{2,1}];
   C2{i}=P2{2,i};          
end
Pop2=(cell2mat(Pop2))';
Cost2=(cell2mat(C2))';

n=numel(pop3);
P3=struct2cell(pop3);
for i=1:n
   PP3{i}=struct2cell(P3{1,i});
   Pop3{i}=[PP3{1,i}{1,1} 
                  PP3{1,i}{2,1}];
    C3{i}=P3{2,i};
end
Pop3=(cell2mat(Pop3))';
Cost3=(cell2mat(C3))';

% xlswrite('pop1-hybrid',Pop1,1);
% xlswrite('pop2-hybrid',Pop2,1);
% xlswrite('pop3-hybrid',Pop3,1);
% xlswrite('pop1-hybrid',Cost1,2);
% xlswrite('pop2-hybrid',Cost2,2);
% xlswrite('pop3-hybrid',Cost3,2);