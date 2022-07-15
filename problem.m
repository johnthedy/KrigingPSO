function [mu,sigma,nRV,dist]=problem(probtype)
if probtype==1
    mu=[0 0];sigma=[1 1];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==2
    mu=[0 0];sigma=[1 1];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==3
    mu=[0 0];sigma=[1 1];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==4
    mu=[0 0];sigma=[1 1];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==5
    mu=[0 0];sigma=[1 1];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==6
    mu=[10 10];sigma=[3 3];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==7
    mu=[0 0];sigma=[1 1];nRV=2;dist={'Normal' 'Normal'};
elseif probtype==8
    mu=[1 1 0.1 0.5 1 0.8];sigma=[0.05 0.1 0.1 0.1 0.15 0.15].*mu;nRV=6;dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
elseif probtype==9
    mu=[1 1 0.1 0.5 1 0.75];sigma=[0.05 0.1 0.1 0.1 0.15 0.15].*mu;nRV=6;dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
elseif probtype==10
    mu=[1 1 0.1 0.5 1 0.6];sigma=[0.05 0.1 0.01 0.05 0.2 0.2];nRV=6;dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
elseif probtype==11
    mu=[200 200 200 50];sigma=mu.*[0.15 0.15 0.15 0.4];nRV=4;dist={'Lognormal' 'Lognormal' 'Lognormal' 'Lognormal'};
elseif probtype==12
    mu=[20000 12 0.04 2e10 9.82e-4 1e11];sigma=mu.*[0.07 0.01 0.12 0.06 0.06 0.06];nRV=6;dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
elseif probtype==13
    mu=[0 0 0 0 0 0 0 0 0 0];sigma=[1 1 1 1 1 1 1 1 1 1];nRV=10;dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal'};
elseif probtype==14
%     mu=[5 42 120 60 3 3 12 90 175 5 10];sigma=[0.1 0.5 1.2 0.6 0.3 0.3 1.2 9 17.5 0.25 0.5];nRV=11;dist={'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Normal' 'Extreme Value' 'Normal' 'Normal' 'Normal' 'Normal'};
    mu=[5 42 119.75 59.75 3 3 12 90 220];sigma=[0.1 0.5 120.25 60.25 0.3 0.3 1.2 9 22];nRV=9;dist={'Normal' 'Normal' 'Uniform' 'Uniform' 'Normal' 'Normal' 'Extreme Value' 'Normal' 'Normal'};
end
end