function x=bound(x,beta1,beta2)
R=rssq(x,2);

for i=1:length(R)
    if R(i)<beta2
        x(i,:)=x(i,:).*(beta2/R(i));
    end
    if R(i)>beta1
        x(i,:)=x(i,:).*(beta1/R(i));
    end
end
end