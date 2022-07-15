function [bestFitness,bestOrganism]=PSO1(ite,particle,x,beta1,beta2,nRV,xdoe,mode,dmin,CP,dmodel)
%% Parameter Setting
%initial velocity
n=nRV;
velocity=ones(particle,n)*0;
weight=0.95;
c1=0.2;
c2=0.2;

%% Start
fitness=fobjPSO(x,xdoe,mode,dmin,CP,dmodel);

[minval,row]=min(fitness);
f_hist(:,1)=fitness;

x_gbest=x(row,:);
minval_history(1)=minval;
x_pbest=x;

for z=2:ite
    %Velocity
    r1=repmat(rand(1,n),particle,1);
    r2=repmat(rand(1,n),particle,1);
    velocity=weight.*velocity+(x_pbest-x).*r1.*c1+(repmat(x_gbest,particle,1)-x).*r2.*c2;
    
    x=x+velocity;
    x=bound(x,beta1,beta2);
    
    for i=1:particle
        if rand()<0.1
            dummy1=unifrnd(-1,1,1,n);
            x(i,:)=dummy1.*(unifrnd(beta2,beta1,1,1)./sqrt(sum(dummy1.^2,2)));
        end
    end
    
    x(1,:)=x_gbest;
    
    %fitness
    fitness=fobjPSO(x,xdoe,mode,dmin,CP,dmodel);
    
    %Record maximum fitness for every generation
    [minval,row]=min(fitness);
    f_hist(:,z)=fitness;
    
    x_gbest=x(row,:);
    minval_history(z)=minval;
    
    for i=1:particle
        if f_hist(i,z)<f_hist(i,(z-1))
            x_pbest(i,:)=x(i,:);
        end
    end
    
end
bestFitness=min(fitness);
bestOrganism=x_gbest;

end