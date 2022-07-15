function [xdoe,G_xdoe,nKs1]=LHS(nRV,nLHSs,nKs1,uxdoe,lxdoe,probtype,mu,sigma,dist)
dummy3=ones(1,nRV).*0.5;
xdoe=lxdoe+dummy3.*(uxdoe-lxdoe);
[sample,~,~,~]=summonsample(1,mu,sigma,nRV,dist,xdoe);
G_xdoe=G(sample,probtype);

stopcrit=0;
while stopcrit<nLHSs
    dummy2=linspace(0,1,nKs1+1);
    [dummy7,~]=size(dummy3);
    
    dummy1=zeros(nKs1-dummy7,nRV);
    for i=1:nRV
        dummy4=randperm(nKs1)';
        dummy6=zeros(dummy7,1);
        for j=1:dummy7
            dummy8=dummy3(j,i)-dummy2;
            [~,dummy5]=min(dummy8(dummy8>=0));
            dummy6(j,:)=dummy5;
        end
        [idx,~]=ismember(dummy4,dummy6);
        dummy4(idx)=[];
        dummy1(:,i)=dummy4;
    end
    
    newdummy3=zeros(nKs1-dummy7,nRV);
    for i=1:nKs1-dummy7
        for j=1:nRV
            newdummy3(i,j)=unifrnd(dummy2(dummy1(i,j)),dummy2(dummy1(i,j)+1));
        end
    end
    
    dummy3=vertcat(dummy3,newdummy3);
    newxdoe=lxdoe+newdummy3.*(uxdoe-lxdoe);
    xdoe=vertcat(xdoe,newxdoe);
    [sample,~,~,~]=summonsample(nKs1-dummy7,mu,sigma,nRV,dist,newxdoe);
    newG_xdoe=zeros(nKs1-dummy7,1);
    for i=1:nKs1-dummy7
        newG_xdoe(i,:)=G(sample(i,:),probtype);
    end
    G_xdoe=vertcat(G_xdoe,newG_xdoe);
    
    stopcrit=length(G_xdoe(G_xdoe<Inf));
    nKs1=nKs1*2;
end
nKs1=nKs1/2;
end
