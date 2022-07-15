if nRV==2 && plotfig==1
    h=figure;
    x0=10;y0=10;width=400;height=380;
    set(gcf,'position',[x0,y0,width,height])

    [~,~,Input1,Input2]=summonsample(1,mu,sigma,nRV,dist,xdoe(1,:));
    [gx1,gx2]=meshgrid(lxdoe-1:0.1:uxdoe+1,lxdoe-1:0.1:uxdoe+1);gx{1}=gx1;gx{2}=gx2;
    dummy0={};
    [a,b]=size(gx1);
    for j=1:nRV
        for k=1:a
            for l=1:b
                dummy0{j}(k,l)=icdf(dist{j},normcdf(gx{j}(k,l),0,1),Input1(j),Input2(j));
            end
        end
    end
    rgx1=dummy0{1};rgx2=dummy0{2};
    if probtype==1
        k=4;m=7;
        z=min(min(0.1.*(rgx1-rgx2).^2-(rgx1+rgx2)./sqrt(2)+k,0.1.*(rgx1-rgx2).^2+(rgx1+rgx2)./sqrt(2)+k),min(rgx1-rgx2+m/sqrt(2),-rgx1+rgx2+m/sqrt(2)));
    elseif probtype==2
        k=5;m=8;
        z=min(min(0.1.*(rgx1-rgx2).^2-(rgx1+rgx2)./sqrt(2)+k,0.1.*(rgx1-rgx2).^2+(rgx1+rgx2)./sqrt(2)+k),min(rgx1-rgx2+m/sqrt(2),-rgx1+rgx2+m/sqrt(2)));
    elseif probtype==3
        k=5;m=9;
        z=min(min(0.1.*(rgx1-rgx2).^2-(rgx1+rgx2)./sqrt(2)+k,0.1.*(rgx1-rgx2).^2+(rgx1+rgx2)./sqrt(2)+k),min(rgx1-rgx2+m/sqrt(2),-rgx1+rgx2+m/sqrt(2)));
    elseif probtype==4
        z=min(min(0.1.*(rgx1-rgx2).^2-(rgx1+rgx2)./sqrt(2)+3,0.1.*(rgx1-rgx2).^2+(rgx1+rgx2)./sqrt(2)+3),min(rgx1-rgx2+3.5*sqrt(2),-rgx1+rgx2+3.5*sqrt(2)));
    elseif probtype==5
        z=0.5.*(rgx1-2).^2-1.5.*(rgx2-5).^3-3;
    elseif probtype==6
        z=2.5-0.2357.*(rgx1-rgx2)+0.00463.*(rgx1+rgx2-20).^4;
    elseif probtype==7
        z=-0.5*(rgx1-rgx2).^2-(rgx1+rgx2)./sqrt(2)+3;
    end
    
    [z1,~]=predictor([gx1(:) gx2(:)],dmodel);
    z1 = reshape(z1, size(gx1));
    hold on
    [~,~]=contour(gx1,gx2,z,[0 0],'LineWidth',0.5,'Color','k');
    [~,~]=contour(gx1,gx2,z1,[0 0],'LineWidth',0.5,'Color','r');
    
    scatter(xdoe(1:end,1),xdoe(1:end,2),10,'filled')
    legend('Limit state function','Kriging prediction','Training samples')
    
    grid on
    xlabel(['{\it U}_{' int2str(1) '}'], 'Fontsize',12);
    ylabel(['{\it U}_{' int2str(2) '}'], 'Fontsize',12);
    set(gca,'fontname','times')
    hold off
end