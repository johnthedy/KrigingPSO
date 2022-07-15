function y=fobjPSO(x,xdoe,mode,dmin,CP,dmodel)
[mu,s2]=predictor(x,dmodel);

w=2;
sig=sqrt(s2);
for i=1:length(mu)
    y(i,:)=(mu(i)*(1-2*normcdf(mu(i)/sig(i)))+sig(i)*(w-sqrt(2/pi)*exp(-0.5*(mu(i)/sig(i))^2)));
end

for i=1:length(mu)
    fp(i,:)=prod(normpdf(x(i,:),0,1));
end

y=1e9-y.*fp;
for i=1:length(y)
    if mu(i)<=0 && mode==1
        w1=10;w2=1;w3=10;
        y(i)=-1/(w1*abs(mu(i))*abs(rssq(x(i,:),2)-rssq(CP,2))+w2*abs(rssq(x(i,:),2)-rssq(CP,2))+w3*length(find(rssq((x(i,:)-xdoe),2)<dmin*3)));
    end
end

for i=1:length(y)
    if min(rssq((x(i,:)-xdoe),2))<dmin
        y(i)=Inf;
    end
end

end

