clear
clc
close all
x=load('day2_0917.txt');
%% ECG signal
y=x(1:95000,1); % ECG signal
t1=linspace(0,95,95000);
figure,plot(t1,y);
title('ECG signal');
xlabel('time');
ylabel('amplitude');

%% PPG signal
z=x(200:95000,2); % PPG signal
t2=linspace(0,94.801,94801);
figure
plot(t2,z,'r');
title('PPG signal');
xlabel('time');
ylabel('amplitude');

%% peak detection of ECG
j=1;
n=length(y);
for i=2:n-1
    if y(i)> y(i-1) && y(i)>= y(i+1) && y(i)> 0.45*max(y)
       maxecgval(j)= y(i);
       maxecgpos(j)=t1(i);
       j=j+1;
     end
end
ecg_peaks=j-1;
maxecg_pos=maxecgpos;
figure,plot(t1,y);
hold on
plot(maxecgpos,maxecgval,'*r');
title('ECG peak');

%% HRV
j=1;
for i=1:ecg_peaks-1
    HRV(j)= maxecg_pos(i+1)-maxecg_pos(i);
    j=j+1;
end 

%% Systolic peak detection of PPG
m=1;
n=length(z);
for i=2:n-1
    if z(i)> z(i-1) && z(i)>= z(i+1) && z(i)> 0.45*max(z)
       sppgval(m)= z(i);
       sppgpos(m)=t2(i);
       m=m+1;
     end
end
sppg_peaks=m-1;
sppg_pos=sppgpos;
sppg_val=sppgval;
figure
plot(t2,z)
hold on
plot(sppgpos,sppgval,'*g');
title('ECG & PPG signal');
legend('ECG signal','PPG signal');

%% Distolic peak dection of PPG
w=1
k=length(z);
for i=2:k-1
    if z(i)>z(i-1) && z(i)>=z(i+1) && z(i)<0.5 && z(i)>=0
        dppgval(w)=z(i);
        dppgpos(w)=t2(i);
        w=w+1;
    end
end
dppg_peaks=m-1;
dppg_pos=dppgpos;
dppg_val=dppgval;
figure,plot(t2,z)
hold on
plot(dppgpos,dppgval,'*g');
title('ECG & PPG signal');

%% Max slope location of PPG
k=length(z);
for i=2:k-1
    slope(i)=(z(i)-z(i-1))./0.001;
end
t=1;
for i=2:length(slope)-1
    if slope(i)>slope(i-1) && slope(i)>slope(i+1) && slope(i)>23
        maxslope(t)=slope(i);
        max_slope_locs(t)=t2(i);
        t=t+1;
    end
end
figure
plot(t2(1:2500),slope(1:2500))
hold on
plot(max_slope_locs,maxslope,'*r')

%% PTTp (distance between ECG peak and PPG peak)
pttp=(sppg_pos-maxecg_pos);
figure,stairs(pttp);
title('PTTp');

%% Regress Model
load final
tb1=table(data((1:35),1),data((1:35),2),bs(1:35,1),'VariableNames',{'pttp','HRV','Systolic_bloodpressure'});
tb2=table(data((36:75),1),data((36:75),2),bs(36:75,1),'VariableNames',{'pttp','HRV','Systolic_bloodpressure'});
tb3=table(data((1:35),1),data((1:35),2),bs(1:35,2),'VariableNames',{'pttp','HRV',' Distolic_bloodpressure'});
tb4=table(data((36:75),1),data((36:75),2),bs(36:75,2),'VariableNames',{'pttp','HRV',' Distolic_bloodpressure'});
md1=fitlm(tb1,'linear','RobustOpts','on') % linear regression (Systolic)

md2=fitlm(tb1,'quadratic','RobustOpts','on')% linear regression (Quadratic) (Systolic)

% md3=fitrsvm(tb1,'bloodpressure') % SVM Regression (Systolic)

md4=fitlm(tb3,'linear','RobustOpts','on') % linear regression (Diastolic)

md5=fitlm(tb3,'quadratic','RobustOpts','on')% linear regression (Quadratic) (Diastolic)
% Testing
y=predict(md1,tb2); % predict data by using linerar regression (Systolic)
y2=predict(md2,tb2); % predict data by using linerar regression(Quadratic) (Systolic)
% y3=predict(md3,tb2); % predict data by using SVM regression (Systolic)
y4=predict(md4,tb4); % predict data by using linerar regression (Distolic)
y5=predict(md5,tb4); % predict data by using linerar regression(Quadratic) (Distolic)
% Training
y6=predict(md1,tb1);
y7=predict(md4,tb3);
%% Estimate  Systolic Blood Pressure 
err1=abs(y-bs(36:75,1)); % linear
% err2=abs(y2-bs(36:75,1)); % SVM
figure,plot(bs(1:35,1))
hold on
plot(y6)
legend('predict data','real data')
xlabel('Sample');ylabel('pressure(mmHg)')
figure,plot(bs(1:35,2))
hold on
plot(y7)
legend('predict data','real data')
xlabel('Sample');ylabel('pressure(mmHg)')
figure,plot(y)
hold on
plot (bs(36:75,1),'r')
legend('predict data','real data')
xlabel('Sample');ylabel('pressure(mmHg)')

figure,plot(y2)
hold on
plot (bs(36:75,1),'r')
legend('predict data','real data')
xlabel('Sample');ylabel('pressure(mmHg)')

%% Estimate  Distolic Blood Pressure 
figure,plot(y4)
hold on
plot (bs(36:75,2),'r')
legend('predict data','real data')
xlabel('Sample');ylabel('pressure(mmHg)')
figure,plot(y5)
hold on
plot (bs(36:75,2),'r')
legend('predict data','real data')
xlabel('Sample');ylabel('pressure(mmHg)')
%% SVM Regression
% md3=fitrsvm(tb1,'bloodpressure')
% y3=predict(md3,tb2) % predict data by using SVM regression
% err3=abs(y3-bs(36:75,1))
% tb3=table(data((1:35),1),data((1:35),2),bs(1:35,2),'VariableNames',{'pttp','HRV','bloodpressure'});
% md4=fitlm(tb3,'linear','RobustOpts','on')
% y4=predict(md4,tb3)