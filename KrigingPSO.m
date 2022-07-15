function [Pf,FE]=KrigingPSO(probtype,uxdoe,lxdoe,beta1,beta2,Kriginginitial,PfCOV,Ulimit,REIFlimit,plotfig)

%% Initial Parameter
[mu,sigma,nRV,dist]=problem(probtype);
modelimit=0.3;
mode=2;
CP=0;expansion=0;levelindex=1;stopcrit=1;
ite=nRV*5;ecosize=nRV*15;dmin=0.01;
nlevel=ceil(nRV*2);
REIFhist=[];Uhist=[];recCP=[];recFELS=[];

%% Initial Kriging using small sample size
[xdoe,G_xdoe,nKs1]=LHS(nRV,1,Kriginginitial,uxdoe,lxdoe,probtype,mu,sigma,dist);
FE=nKs1;
regfunc=@regpoly0; corrfunc=@corrgauss;
theta = ones(1,size(dist,1));
lob = ones(1,size(dist,1))*1e-2; upb = ones(1,size(dist,1))*1e2;
dmodel=dacefit(xdoe,G_xdoe,regfunc,corrfunc,theta,lob,upb);

%% Looping to strengthen Kriging
indexloop=1;
while true
    %Perform line search technique
    if isempty(G_xdoe(G_xdoe<0))~=1 && expansion==0
        [beta1,beta2,recCP,CP,FELS]=linesearch2(xdoe,G_xdoe,recCP,expansion,levelindex,nlevel,dmodel);
        recFELS(indexloop)=FELS;
    elseif isempty(G_xdoe(G_xdoe<0))~=1 && expansion==1
        if mod(indexloop,5)==0
            [beta1,beta2,recCP,CP,FELS]=linesearch2(xdoe,G_xdoe,recCP,expansion,levelindex,nlevel,dmodel);
            recFELS(indexloop)=FELS;
        end
    end
    
    %Check to enter leveling stage
    error=1.5;nDATA=4;
    if expansion==0 && length(G_xdoe(G_xdoe<0))>nRV*2 && length(recCP)>nDATA && 100*(abs(mean(recCP(end-nDATA:end))-recCP(end)))/recCP(end)<error && 100*(abs(max(recCP(end-nDATA:end))-recCP(end)))/recCP(end)<error && 100*(abs(min(recCP(end-nDATA:end))-recCP(end)))/recCP(end)<error
        expansion=1;modelimit=0.85;finalCP=recCP(end);
    end
    
    %Check if require go back to CP stage
    if expansion==1 && 100*abs(finalCP-CP)/finalCP>10
        expansion=0;modelimit=0.3;beta1=CP+3;beta2=CP;stopcrit=0;
    end
    
    %Check criteria to enter next level
    nDATA=max(10,nRV^2-15);
    if stopcrit>1 || (length(Uhist)>nDATA && min(Uhist(end-nDATA:end))>5 && max(REIFhist(end-nDATA:end))<0 && isempty(G_xdoe(G_xdoe<0))==0)
        levelindex=levelindex+1;stopcrit=0;
    end
    
    %Determine objective function type that will be used in PSO
    if rand()<modelimit || length(G_xdoe(G_xdoe<0))<1
        mode=2;
    else
        mode=1;
    end
    
    %Perform PSO to find best particle
    dummy1=unifrnd(-1,1,ecosize,nRV);
    eco=dummy1.*(unifrnd(beta2,beta1,ecosize,1)./sqrt(sum(dummy1.^2,2)));
    [~,bestOrganism]=PSO1(ite,ecosize,eco,beta1,beta2,nRV,xdoe,mode,dmin,CP,dmodel);
    
    %Calculate and record best particle U function and REIF function
    %Uhist is matrix that record U function value for each obtained best particle
    %REIFhist is matrix that record REIF function value for each obtained best particle
    [y,s2]=predictor(repmat(bestOrganism,2,1),dmodel);
    y=y(1);sig=sqrt(s2(1));w=2;
    if mode==2
        Uhist=vertcat(Uhist,abs(y)/sig);
        REIFhist=vertcat(REIFhist,(y*(1-2*normcdf(y/sig))+sig*(w-sqrt(2/pi)*exp(-0.5*(y/sig)^2))));
        if (Uhist(end)>Ulimit && expansion==1) || (REIFhist(end)<=REIFlimit && expansion==1)
            stopcrit=stopcrit+1;
        end
    end
    
    %Update Kriging doe
    xdoe=vertcat(xdoe,bestOrganism);
    [sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,bestOrganism);G_xdoe=vertcat(G_xdoe,G(sample,probtype));FE=FE+1;
    dmodel=dacefit(xdoe,G_xdoe,regfunc,corrfunc,theta,lob,upb);
    
    %Check stopping criteria
    if levelindex==nlevel+1
        break
    end
    indexloop=indexloop+1;
end

%% Calculate Pf
sampleincrement=1e5;index=0;rsamplesize=Inf;totalevalsample=0;
while sampleincrement*index<rsamplesize
    totalevalsample=totalevalsample+sampleincrement;
    CDFthreshold=chi2cdf((0.9*CP)^2,nRV);
    CDF=unifrnd(CDFthreshold,1,sampleincrement,1);
    sample=normrnd(0,1,sampleincrement,nRV);R_xxx=rssq(sample,2);R_chi2=sqrt(chi2inv(CDF,nRV));
    sample=sample.*R_chi2./R_xxx;
    [y,~]=predictor(sample,dmodel);
    index=index+1;
    temp(index)=length(y(y<0));
    Pf=sum(temp)/(sampleincrement*index)*(1-CDFthreshold);
    rsamplesize=(1-CDFthreshold)*(1-Pf)/(PfCOV*PfCOV*Pf);
end
time1=toc;
fprintf( 'Computation time = %s second\n',time1)
clear temp

%% Plot figure
plotfigure
time2=toc;
fprintf( 'Plotting time = %s second\n',time2-time1)

%% Print Result
fprintf( 'Failure Probability Result = %d\n',Pf)
fprintf( 'Total Real Function Evaluation = %d\n',FE)
fprintf( 'Total Kriging Function Evaluation during Kriging loop = %d\n',sum(recFELS))
fprintf( 'Total Kriging Function Evaluation during Pf evaluation = %d\n',totalevalsample)
end
