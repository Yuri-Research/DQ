function [dqF,isomodel]=performDQmv(r,y,u,display,specifications,values)
%This function performs data quality assessment using
%The specifications represent the named variables that we wish to change
%(those whose values are given as global)
%and Values represents the desired value. Both Specifications and Values
%should be a cell arrays.
%Variable names of interest:
%alpha1, Nis, Nd, thetaD, lambdamy, lambdany, lambdamL, lambdanL,...
%lambdamu, lambdanu, etaL1, etaY, etaCond, lambdaR, etaChi, lambdaRy,...
%etaCondy, etaChiy, etaUy
%This is the multivariate input version of the programme previously written

%Copyright 2015, 2022 Yuri A.W. Shardt

global alpha1 Nis; %Laguerre model parameters
global Nd thetaD; %ARX based parameters
global lambdamy lambdany lambdamL lambdanL lambdamu lambdanu etaL1 etaY etaCond lambdaR etaChi; %forgetting values and limits for Isaakson's approach
global display2;
display2=display;

%Default free parameter values
alpha1=0.8; %time constant of Laguerre model
Nis=6; %Order of the Laguerre model

Nd=5; %Number of orders in the ARX model
thetaD=20; %Time delay (assumed)

%Model forgetting values and limits for the Laguerre-based approach
lambdamy=0.99;
lambdany=0.99;
lambdamL=0.99;
lambdanL=0.99;
lambdamu=0.99;
lambdanu=0.99;
etaL1=1e-7;
etaY=1e-10;
etaCond=1e6;
lambdaR=0.95;
etaChi=10^3;

%Overwrite any specified parameter values
if exist('specifications','var')
    a=length(specifications);
    for i=1:a
        eval([specifications{i} '=' values{i} ';']);
    end
end

dqF=NaN;
%Check that the data makes sense--dimension and lengths:
nR=size(r);
nY=size(y);
nU=size(u);

if length(nR) > 2 || length(nY) > 2 || length(nU) > 2
    fprintf('At least one of the three matrices entered have a dimension greater than 2.\n');
    return;
end

nR_1=max(nR);
nR_2=min(nR);

nY_1=max(nY);
nY_2=min(nY);

nU_1=max(nU);
nU_2=min(nU);

if nR_1 ~= nY_1 || nR_1 ~= nU_1 || nY_1 ~= nU_1
    fprintf('At least one of the three matrices entered are not the same length as the others.\n');
    return;
end

if nR_2 ~= 1 || nY_2 ~= 1
    fprintf('At least one of the reference or output signals are not univariate.\nUnfortunately, the programmes as it currently stands can only take univariate cases for these two signals.\n');
    return;
end

%Scale the data
rScaled=(r-nanmean(y))./nanstd(y); %The scaling on the setpoint and output should be the same, so I am using the output to determine the scaling for both
yScaled=(y-nanmean(y))./nanstd(y);
if (nanstd(u)~=0)
    uScaled=u;%zscore(u);
else
    uScaled=u;
end

%Determine the modes
k_mode=1;
while k_mode<=length(u)
    if k_mode==1
        mode(k_mode)=determineMode(y(k_mode),NaN,r(k_mode+1));        
    elseif k_mode==length(r)
        mode(k_mode)=determineMode(y(k_mode),r(k_mode-1),NaN);
    else
        mode(k_mode)=determineMode(y(k_mode),r(k_mode-1),r(k_mode+1));
    end
    if isnan(u(k_mode))
            mode(k_mode)=-1;
    end
    k_mode=k_mode+1;
end

%for each mode compute the appropriate test metric
locatechange=find(diff(mode)~=0);
qLast=1;
if isempty(locatechange)
    %Only a single region was found.
    display2.Message='Processing Single Mode';
    [partition,isomodel]=selectMethod(rScaled,yScaled,uScaled,mode(qLast));
else
    qCurrent=mode(locatechange(1));

    partition=-0*ones(length(u),1); %determine the partitions
    isomodel=NaN*ones(length(u),1); %determine whether the model is the same--only applies for closed-loop without r_t for now.

    for regions=1:length(locatechange)
        qCurrent=locatechange(regions);
       
        display2.Message=['Processing Mode ' int2str(regions) ' of ' int2str(length(locatechange)) ' modes.'];
        if mode(qCurrent)==-1
            %We have bad data.
            if qLast>1
                partition(qLast:qCurrent)=0+partition(qLast-1);
            else
                partition(qLast:qCurrent)=0;
            end
            
        else
    
        if qCurrent-qLast>20;
            %We need at least a good sized data set, assume that 20 samples is
            %the bare minimum.
            %Select the method based on the mode.
            if qLast>2
                [partition(qLast:qCurrent),isomodel(qLast:qCurrent)]=selectMethod(uScaled(qLast:qCurrent,:),yScaled(qLast:qCurrent),uScaled(qLast:qCurrent,:),mode(qLast));
                partition(qLast:qCurrent)=partition(qLast:qCurrent)+partition(qLast-1);

            else
                [partition(qLast:qCurrent),isomodel(qLast:qCurrent)]=selectMethod(uScaled(qLast:qCurrent,:),yScaled(qLast:qCurrent),uScaled(qLast:qCurrent,:),mode(qLast));
            end
            if display2.CancelRequested
                fprintf('Cancelled the computation.\n')
                partition=0;
                return;
            end

        end
        end
        qLast=qCurrent+1;
    end
    if qLast<length(u)
        %One last run to deal with the final list needs to be done
        qLast=qCurrent+1;
        qCurrent=length(u);
        if qCurrent-qLast>20;
            %We need at least a good sized data set, assume that 20 samples is
            %the bare minimum.
            %Isaakson's Method
            [partition(qLast:end),isomodel(qLast:end)]=selectMethod(rScaled(qLast:qCurrent),yScaled(qLast:qCurrent),uScaled(qLast:qCurrent,:),mode(qLast));
            partition(qLast:end)=partition(qLast:end)+partition(qLast-1);
        end
    end
end

Nent=20;

dqF=simplifyDataPartitions(partition);

return
%% Method to select the correct method given the operating mode
function [partitioned,entropy2]=selectMethod(r,y,u,mode)
%Select and partition the data based on the desired method. Entropy2 is
%only returned if the data is closed-loop with no changes in r_t.
global Nd thetaD;
entropy2=-inf*ones(length(u),1);
switch mode
    case 0
        partitioned=isaaksons(r,y,u); %Open-loop
        partitioned=reduce(partitioned,y,u);
    case -1
        %Do nothing with the bad data
        partitioned=-500*ones(size(y));
    otherwise
        partitioned=-500*ones(size(y));
end

return
%% Reduce the number of partitions
function reducedP=reduce(partitioning,y,u)
changeD=abs(diff(diff(partitioning)));
locations=find(changeD);
nLoc=length(locations);
count=1;
curEnt=0;
prevEnt=NaN;
reducedP=partitioning;
entropyVal=[];
[u1,u2]=size(u);
while count<nLoc
    start1=locations(count);
    end1=locations(count+1);
    if abs(end1-start1)<5
        %Ignore all regions less than 10 in size
        count=count+2;
        continue;
    end
    entrY((count-1)/2+1)=entropy(y(start1:end1));
    entrU((count-1)/2+1)=0;
    for i=1:u2
    entrU((count-1)/2+1)=entrU((count-1)/2+1)+entropy(u(start1:end1,i));
    end
    curEnt=entropy(y(start1:end1))-entropy(u(start1:end1));
    
    if abs(curEnt-prevEnt)<0.25
        %consider regions to be the same
        %merge them.
        if count>1
            diff2=partitioning(start1)-partitioning(locations(count-1));
        else
            diff2=0;
        end
        reducedP(locations(count-1)+1:end1)=reducedP(locations(count-1)); %Take the previous count number and use it. We should now decrease the remaining values by the difference
        reducedP(end1+1:end)=reducedP(end1+1:end)-diff2;
        
        prevEnt=(prevEnt+curEnt)/2; %Obtain a better estimate of the entropy for the region        
    else
        prevEnt=curEnt;
    end
    entropyVal((count-1)/2+1)=curEnt;
    count=count+2;
end

return

%% Determine the mode given the data
function mode=determineMode(y,a,b)
%Based on the value of rt determine the current mode.
mode=0; %normal operation
if abs(y)<0.0001 %To avoid certain issues, I am making values close to zero be equal to zero and hence process upsets.
    mode=-1; %bad data
elseif isnan(y)
    mode=-1;
else
    %There are currently no other choices.
end
return

%% Implement Isaakson's method for data segmentation
function partition=isaaksons(rScaled,yScaled,uScaled)
partition=-1*ones(length(uScaled),1); %This variable keeps track of the current parition and to which region each data set should be assigned.
partitioncount=1; %This variable counts what the current partition count is.
global alpha1 Nis;
global lambdamy lambdany lambdamL lambdanL lambdamu lambdanu etaL1 etaY etaCond lambdaR etaChi
global display2

k=1;
while (sum(uScaled(k,:)) == 0 && k<length(uScaled))
    k=k+1;
end

%Check if a value was found
if k>=length(uScaled)
    fprintf('No non-zero inputs found. Procedure failed\n');
    return;
end
%Compute parameters required for Isaakson's method
%Create Laguerre filters
%Using global to pass the values, so that I can change them dynamically.
for i=1:Nis
    Lmodel{i}=tf([sqrt(1-alpha1)],[1 -alpha1],1)*(tf([-alpha1 1],[1 -alpha1],1))^(i-1);
end
diff=-0.35;
for i=1:Nis
    Lmodel2{i}=tf([sqrt(1-alpha1-diff)],[1 -alpha1-diff],1)*(tf([-alpha1-diff 1],[1 -alpha1-diff],1))^(i-1);
end
[u1,u2]=size(uScaled);
counter=1;
for i=1:Nis
    
    for j=1:u2
        if j==1
    L(:,counter)=lsim(Lmodel{i},uScaled(:,u2));
        else
            L(:,counter)=lsim(Lmodel2{i},uScaled(:,u2));
        end
    counter=counter+1;
    end
end

my=0;
mL1=0;
mU=0;



nuY(k)=0;
nuL1(k)=0;
nuU(k)=0;

mL1_all=zeros(1,u2);
mU_all=zeros(1,u2);

for i=k+1:length(uScaled)
    my=lambdamy*yScaled(i)+(1-lambdamy)*my;
    nuY(i)=(2-lambdamy)/2*(lambdany*(yScaled(i)-my))^2+(1-lambdany)*nuY(i-1);
    nuL1(i)=0;
    nuU(i)=0;
    for l=1:u2
        mL1=mL1_all(l);
        mU=mU_all(l);
        mL1_all(l)=lambdamL*L(i,1)+(1-lambdamL)*mL1;
        nuL1(i)=nuL1(i)+(2-lambdamL)/2*(lambdanL*(L(i,1)-mL1))^2+(1-lambdanL)*nuL1(i-1);
        mU_all(l)=lambdamu*uScaled(i)+(1-lambdamu)*mU;
        nuU(i)=nuU(i)+(2-lambdamu)/2*(lambdanu*(uScaled(i)-mU))^2+(1-lambdanu)*nuU(i-1);
    end

    
end
kest=-1;
%R{={};

for i=k:length(uScaled) 
    if kest==-1
        kest=i;
        R{i}=zeros(1,1);
    else
        R{i}=lambdaR*R{i-1}+(1-lambdaR)*L(i,:)*L(i,:)';
    end


    PhiK=[L(kest:i,:) yScaled(kest:i)];
    [q,r]=qr(PhiK);
    if i-kest-Nis<=0
        chiN(i)=NaN;
        partitioncount=partitioncount+1;
        partition(i)=partitioncount;
        continue;
    end
    R1=r(1:Nis,1:Nis);
    R2=r(1:Nis,Nis+1);
    R3=r(Nis+1,Nis+1);
    chiN(i)=norm(R2*sqrt(i-kest+1)/abs(R3),2);

    if cond(R{i}) > etaCond
        kest=-1;
        partitioncount=partitioncount+1;
        partition(i)=partitioncount;
        continue;
    end
    if chiN(i)>etaChi
        kest=-1;
        partitioncount=partitioncount+1;
        partition(i)=partitioncount;
        continue;
    end
display2.Value=i/length(uScaled);

if display2.CancelRequested
    partition=NaN;
    fprintf('User stopped the programme\n');
    return
end
    partition(i)=partitioncount;
end
return
%% Entropy-based spurious partition removal
function [slope,r2]=entropy(signal)
L1=zeros(length(signal)-1,1);
nW=length(signal);
for i=2:length(signal)
    if i==2
        L(2)=abs(signal(i-1)-signal(i));
    else
        L(i)=abs(signal(i-1)-signal(i))+L(i-1);
        if i>nW
            L1(i)=L(i)-L(i-nW);
        else
            L1(i)=L(i);
        end
    end
end
[a,b,c,d,e]=regress(L1,[ones(length(L1),1) [0:1:length(L1)-1]']);
slope=a(2);
r2=e(1);


