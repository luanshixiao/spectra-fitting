clear
close all

%% read folder

% path
path="ADD YOUR FOLDER PATH HERE"+"\";
files=dir(path+"*.asc");
% film thickness(nm)
d=700;

%% read opts

opts = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines = [91, Inf];
opts.Delimiter = "\t";
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "skip";

%% for

A=zeros(500,6,length(files));
for ii=1:length(files)
    
    %% prepare spectrum data
    
    filename=files(ii).name;
    sample(ii)=string(filename(1:end-4));
    spec=readtable(path+filename,opts);
    spec=table2array(spec);
    
    T_all=spec(:,2)/100;
    % exclude T<0
    T_all(T_all<0.01)=0;
    % exclude T>0
    T_all(T_all>1)=1;
    
    %% convert data
    
    % lambda
    lambda=spec(:,1);
    % x_all=E=hv
    x_all=1240./lambda;
    % T_all=exp(-alpha_all*d) >> alpha_all=-ln(T_all)/d
    alpha_all=-log(T_all)./(d*1e-7);
    % y_all=sqrt(alpha_all*hv)
    y_all=sqrt(alpha_all.*x_all);
    
    %% fitting region
    
    % lower
    first=find(alpha_all>20000);
    first=first(1);
    % upper
    last=find(alpha_all<50000);
    last=last(end);
    
    x=x_all(first:last);
    y=y_all(first:last);
    
    
    %% fit
    
    x0=x(1):0.001:x(end);
    y0=interp1(x,y,x0,'linear');
    
    [xData, yData] = prepareCurveData( x0, y0 );
    
    region=find(xData>(xData(1)+0.15));
    region=region(1);
    step=floor((region-1)/10);
    
    jj=1;
    while jj>=1
        
        left=1+(jj-1)*step;
        right=left+region;
        Eleft=xData(left);
        Eright=xData(right);
        xfit=xData(left:right);
        yfit=yData(left:right);
        ft = fittype( 'poly1' );
        [fitresult, gof] = fit( xfit, yfit, ft );
        
        A(jj,1,ii)=left;
        A(jj,2,ii)=right;
        A(jj,3,ii)=Eleft;
        A(jj,4,ii)=Eright;
        A(jj,5,ii)=gof.adjrsquare;
        A(jj,6,ii)=gof.rmse;
        
        % fitting linestyle: y=kx+b
        eq(jj,1)=fitresult.p1;%k=p1
        eq(jj,2)=fitresult.p2;%b=p2
        
        if (right+step)>length(xData)
            break
        end
        
        jj=jj+1;
        
    end
    
    %% the best fitresult
    
    [~,r]=max(A(:,5,ii));
    k=eq(r,1);
    b=eq(r,2);
    
    %% calculate Eg and slope by fitresult
    
    Eg(ii)=-b/k;
    slope(ii)=k;
    
    %% plot
    
    figure('Visible','off')
    hold on
    
    f1=plot(x_all,y_all);
    f1.LineWidth=2;
    f1.Color='black';
    
    f2=fplot(@(x) k*x+b);
    f2.LineWidth=1;
    f2.Color='blue';
    f2.LineStyle='--';
    
    t1=text(3,100,"Eg="+num2str(roundn(Eg(ii),-2)));
    t1.FontSize=20;
    
    t2=text(3,50,"slope="+num2str(roundn(slope(ii),0)));
    t2.FontSize=20;
    
    axis([1.6 3.6 0 500])
    xlabel('$$ h{\nu} $$','Interpreter','latex');
    ylabel('$$ ({\alpha}h{\nu})^{1/2} $$','Interpreter','latex');
    
    title(sample(ii))
    
    %% save plot
    
    
    print(path+sample(ii),'-dpng')
end

%% result table

ret=table(Eg',slope','VariableNames',{'Eg','slope'},'RowNames',sample');
writetable(ret,path+"Eg&slope.csv",'Delimiter',',','WriteRowNames',true)
