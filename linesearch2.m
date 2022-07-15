function [beta1,beta2,recCP,CP,FELS]=linesearch2(xdoe,G_xdoe,recCP,expansion,levelindex,nlevel,dmodel)
dummy8=[];FELS=0;
Rlist=rssq(xdoe,2);
for j=2:length(G_xdoe)
    dummy1=Rlist(j);
    if G_xdoe(j)<0
        dummy3=G_xdoe(j);
        while dummy3<0
            dummy1=dummy1-0.1;
            dummy2=xdoe(j,:).*dummy1/rssq(xdoe(j,:),2);
            [y,~]=predictor(repmat(dummy2,2,1),dmodel);
            FELS=FELS+1;
            dummy3=y(1);
        end
        dummy8=vertcat(dummy8,dummy1);
    else
        dummy3=G_xdoe(j);
        dummy2=xdoe(j,:).*7/rssq(xdoe(j,:),2);
        [y,~]=predictor(repmat(dummy2,2,1),dmodel);
        if y(1)<0 && dummy1<=7
            while dummy3>0
                dummy1=dummy1+0.1;
                dummy2=xdoe(j,:).*dummy1/rssq(xdoe(j,:),2);
                [y,~]=predictor(repmat(dummy2,2,1),dmodel);
                FELS=FELS+1;
                dummy3=y(1);
                if dummy1>7
                    break
                end
            end
            dummy8=vertcat(dummy8,dummy1);
        end
    end
end

[CP,~]=min(dummy8);
level=linspace(CP*0.95,CP+3,nlevel+1);
recCP=vertcat(recCP,CP);
if isempty(recCP)~=1
    beta1=CP+3;beta2=CP*0.95;
end
if expansion~=0
    beta1=level(levelindex+1);
    beta2=CP*0.95;
end
end