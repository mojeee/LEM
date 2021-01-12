% Master project 05/19/2019
% clear all;
clc;
cd 'F:\test\DNS_cases\LaminarFlameSpeedNewMechanism\LEMNewFrequencyEquation82(copy)\Turbulence_K74';
% Number of Profile saved
NoPs=12450;

% read all data of ember code from h5py file
for i=1:NoPs
    
    if i<=9
        string1 = "prof00000";
        number = num2str(i);
        string2 = ".h5";
        profile=string1+number;
        profile2=profile+string2;
        Temperature(:,i)=h5read(profile2,'/T');
%         species(:,:,i)=h5read(profile2,'/Y');
        position=h5read(profile2,'/x');
    
    
    elseif i>9 && i<=99
            
        string1 = "prof0000";
        number = num2str(i);
        string2 = ".h5";
        profile=string1+number;
        profile2=profile+string2;
        Temperature(:,i)=h5read(profile2,'/T');
%         species(:,:,i)=h5read(profile2,'/Y');
        position=h5read(profile2,'/x');
        
     elseif i>99 && i<=999
         
         string1 = "prof000";
         number = num2str(i);
         string2 = ".h5";
         profile=string1+number;
         profile2=profile+string2;
         Temperature(:,i)=h5read(profile2,'/T');
%          species(:,:,i)=h5read(profile2,'/Y');
         position=h5read(profile2,'/x');
     
    elseif i>999 && i<=9999
         
         string1 = "prof00";
         number = num2str(i);
         string2 = ".h5";
         profile=string1+number;
         profile2=profile+string2;
         Temperature(:,i)=h5read(profile2,'/T');
%          species(:,:,i)=h5read(profile2,'/Y');
         position=h5read(profile2,'/x');
     elseif i>9999 && i<NoPs
         
         string1 = "prof0";
         number = num2str(i);
         string2 = ".h5";
         profile=string1+number;
         profile2=profile+string2;
         Temperature(:,i)=h5read(profile2,'/T');
%          species2(:,:,i-9999)=h5read(profile2,'/Y');
         position=h5read(profile2,'/x');
        
    end
    i
end

% Progress variabe calculation 
%T_equilibrium=2000; % K
%T_initial=300; % K

% % Progress Variable
% % for j=1:1127
% %     for i=1:length(Temperature)
% %         PV(i,j)=(Temperature(i,j)-T_initial)/(T_equilibrium-T_initial);
% %     end
% %     
% % end
% 
% % write Progress variable in a csv file
% % using LEM data to create PDFccccccccc
% fid = fopen('interpedc.csv','w');
% for j=1:1125
%     for i=1:length(Temperature)
%         fprintf(fid,'%f, ',Temperature(i,j));
%         fprintf(fid,'\n'); 
%     end
%    
% end
% fclose(fid);

%%

cd 'F:\Thesis\DNS data\K6_k1';

load 'Yi[2].mat';


%%

% density calculation
pressure=101325;
R=8.314;
W=[2.01 1.00 15.99 31.99 17.00 18.01 33.00 34.01 15.03 16.04 28.01 44.00 29.01 30.02 31.03 28.01]; 
W=W/1000;
 meanweight=zeros(3850,19998);
 %%
for i=1:3850
    
    for j=1:19998
        i
        j
        MM=0;
        
        for k=1:16
            
            if j<=9999
            MM= MM+species(i,k,j)/W(1,k);
            elseif j>9999
            MM= MM+species2(i,k,j-9999)/W(1,k);
            end
        end
        meanweight(i,j)=1/MM;
        density(i,j)=(pressure*meanweight(i,j))/(R*Temperature(i,j));

    end
    
    
end
%%
for i=1:3850
    
    for j=1:19998
        
       density(i,j)=(pressure*meanweight(i,j))/(R*Temperature(i,j));
        

    end
    
    
end


%%

cd 'F:\test\Data\A1';

save('RhoLEMlambda2ka6.mat','density');



%%



% so wdotCO2 is 3850 by 20000

for i=1:3850
    
    for j=1:19998
        i
        j
     wdotCO2(i,j)=44.0095*(OHLEM(i,j)*Rho(i,j)/17.00734)*(COLEM(i,j)*Rho(i,j)/28.01010)*1.51*10^(4)*Temperature(i,j)^1.3*exp(758/1.98720426/Temperature(i,j)); 
     wdotCO2(i,j)=wdotCO2(i,j)-44.0095*(HLEM(i,j)*Rho(i,j)/1.007940)*(CO2LEM(i,j)*Rho(i,j)/44.00950)*1.57*10^(6)*Temperature(i,j)^(1.3)*exp(-22337/1.98720426/Temperature(i,j)); 

    end
    
    
end

%%

sum=zeros(12500);
for j=1:12500
    
    for i=1:1032

        sum(j)=sum(j)+wdotCO2(i,j);


    end
    
end

plot (sum);

%%
% mean temperature profile calculation

Tmean=zeros(3850,1);
for i=1:3850
    
    for j=1:19999
        
        Tmean(i,1)=Tmean(i,1)+Temperature(i,j);
        
    end
    
    
    
end
%      p=plot(position(:,1),Temperature(:,1),'k','LineWidth',2);
hold on;
plot(position,Tmean/20000,'k','LineWidth',2);

xlabel('Position (m)');
ylabel('Temperature (k)');







%%

% species1=zeros(3850,6,20000);

for i=1:3850
    for j=1:9999
        HLEM(i,j)=species(i,2,j);

%         H2LEM(i,j)=species(i,1,j);
%         O2LEM(i,j)=species(i,4,j);
%         OHLEM(i,j)=species(i,5,j);
%         CH4LEM(i,j)=species(i,10,j);
%         COLEM(i,j)=species(i,11,j);
%         CO2LEM(i,j)=species(i,12,j);

        
    end
end
    
%     clear species;
    
for i=1:3850
    for j=1:9999
                HLEM(i,j+9999)=species2(i,2,j);

% %         H2LEM(i,j+9999)=species2(i,1,j);
% %         O2LEM(i,j+9999)=species2(i,4,j);
% %         OHLEM(i,j+9999)=species2(i,5,j);
% %         CH4LEM(i,j+9999)=species2(i,10,j);
% %         COLEM(i,j+9999)=species2(i,11,j);
%         CO2LEM(i,j+9999)=species2(i,12,j);
% 
%         
    end
end
    
    
%     clear species2;




%%

for i=1:1032
    
    for j=1:12500
        
        CO2LEM15(i,j)=species(i,12,j);
    end
%     for j=10001:20000
%         
%         CO2LEM1(i,j)=species2(i,12,j-10000);
%     end
end

%%
cd 'F:\test\Data\A2';

save('PosLEMlambda2ka6.mat','position');
save('TempLEMlambda2ka6.mat','Temperature');
save('CH4LEMlambda2ka6.mat','CH4LEM');
save('CO2LEMlambda2ka6.mat','CO2LEM');
save('COLEMlambda2ka6.mat','COLEM');
save('H2LEMlambda2ka6.mat','H2LEM');
save('O2LEMlambda2ka6.mat','O2LEM');
save('OHLEMlambda2ka6.mat','OHLEM');



%%

a=10;
for i=240:399
    
    end1=3850-a;
    for j=1:3850
        Temperature1(j,i)=Temperature(j,i);
    end
    for j=1:end1
        Temperature(j+a,i)=Temperature1(j,i);
    end
    for j=1:3850
        for k=1:16
            species1(j,k,i)=species(j,k,i);
        end
    end
    for j=1:end1
        for k=1:16
            species(j+a,k,i)=species1(j,k,i);
        end
    end
    a=a+3;
    i
end
a=470;
for i=400:999
    
    end1=3850-a;
    for j=1:3850
        Temperature2(j,i)=Temperature(j,i);
    end
    for j=1:end1
        Temperature(j+a,i)=Temperature2(j,i);
    end
    for j=1:3850
        for k=1:16
            species2(j,k,i)=species(j,k,i);
        end
    end
    for j=1:end1
        for k=1:16
            species(j+a,k,i)=species2(j,k,i);
        end
    end
    i
end

%%
clear species2;
clear species1;
clear Temperature2;
clear Temperature1;
clear i;
clear j;
clear k;
clear a;
clear end1;

%%
%initial condition
cd 'F:\test\DNS_cases\LaminarFlameSpeedNewMechanism\LEMNewFrequencyEquation84(copy)';

profile2="prof000062.h5";
Temperatureinitial(:)=h5read(profile2,'/T');
speciesinitial(:,:)=h5read(profile2,'/Y');
positioninitial=h5read(profile2,'/x');
hold on;
plot(positioninitial,Temperatureinitial,'b','LineWidth',2);
xlabel('Position (m)');
ylabel('Temperature (K)');


%%
% subplot(1,2,1)
% for i=500:12500
    hold on;
%         p=plot(position(:,1),Temperature(:,i),'--r');
%         p=plot(position(:,1),Temperature(:,522),'-r','linewidth',1.3);
%         p=plot(position(:,1),Temperature(:,532),'--b','linewidth',1.3);
%         p=plot(position(:,1),Temperature(:,722),'-.k','linewidth',1.3);
%         p=plot(position(:,1),Temperature(:,1619),'-r','linewidth',1.3);
%         p=plot(position(:,1),Temperature(:,1997),'-b','linewidth',1.3);
%         p=plot(position(:,1),Temperature(:,1733),'-m','linewidth',1.3);
%         pause(0.5);
%         delete(p);
%  i
% end

Tmean=zeros(1032,1);
for i=1:1032
    
    for j=4000:12500
        
        Tmean(i,1)=Tmean(i,1)+Temperature(i,j);
        
    end
    
    
    
end
%      p=plot(position(:,1),Temperature(:,1),'k','LineWidth',2);
hold on;
box on;
% title ( 'Case A1');
plot(position,Tmean/8500,'r','LineWidth',1.2);

xlabel('Domain (m)');
ylabel('Temperature (k)');
%%
Tmean=Tmean./7000;
maxT=max(Tmean);
%%
b=1;
for i=1:3850
    
    if abs(Tmean(i,1)-(0.95)*maxT)< 0.1
        point=i;
        break;
    end
    if abs(Tmean(i,1)- 315)< 5 && b==1
        point2=i;
        b=5
    end
end



%%
% major Scatter 


hold on;
subplot(3,2,2);
plot(T,CH4);
xlabel('Temperature (K)');
 ylabel('Y_{CH4}');
subplot(3,2,4);
plot(T,CO2);
xlabel('Temperature (K)');
 ylabel('Y_{CO2}');
subplot(3,2,6);
plot(T,O2);
xlabel('Temperature (K)');
 ylabel('Y_{O2}');
% subplot(2,1,1);

for i = 500:20:12500
%      p=plot(position(:,1),Temperature(:,1),'k','LineWidth',2);
     hold on;
%      p=plot(position(:,1),Temperature(:,i),'--r');
% p=plot(species(:,12,1),species(:,11,1),'k','LineWidth',1.5);
% if i<=10000
    subplot(3,2,1);
%      p=plot(Temperature(:,1),species(:,10,1),'r','LineWidth',1.5);
     hold on;
     plot(Temperature(:,i),CH4LEM(:,i),'k');
     xlabel('Temperature (K)');
 ylabel('Y_{CH4}');
         subplot(3,2,3);
%      p=plot(Temperature(:,1),species(:,12,1),'r','LineWidth',1.5);
     hold on;
     plot(Temperature(:,i),CO2LEM(:,i),'k');
     xlabel('Temperature (K)');
 ylabel('Y_{CO2}');
         subplot(3,2,5);
%      p=plot(Temperature(:,1),species(:,4,1),'r','LineWidth',1.5);
     hold on;
     plot(Temperature(:,i),O2LEM(:,i),'k');
     xlabel('Temperature (K)');
 ylabel('Y_{O2}');
%      p=plot(species(:,12,i),species(:,11,i),'--r');
i
%     
% elseif i>10000
%     subplot(3,2,5);
% %      p=plot(Temperature(:,1),species(:,1,1),'r','LineWidth',1.5);
%      hold on;
%      plot(Temperature(:,i),species2(:,4,i-10000),'k');
% %      xlabel('Temperature (K)');
% %  ylabel('Y_{CO2}');
%          subplot(3,2,3);
% %      p=plot(Temperature(:,1),species(:,11,1),'r','LineWidth',1.5);
%      hold on;
%      plot(Temperature(:,i),species2(:,12,i-10000),'k');
% %      xlabel('Temperature (K)');
% %  ylabel('Y_{CO}');
%          subplot(3,2,1);
% %      p=plot(Temperature(:,1),species(:,5,1),'r','LineWidth',1.5);
%      hold on;
%      plot(Temperature(:,i),species2(:,10,i-10000),'k');
%      xlabel('Y_{CO2}');
%      ylabel('Y_{CO}');
%      title('LEM Results K74');
%      grid on;
     hold on;
%       pause(0.09);
%       delete(p);
% i
% end
 end
    subplot(3,2,1);
     p=plot(Temperature(:,1),CH4LEM(:,1),'r','LineWidth',1.5);
              subplot(3,2,3);
     p=plot(Temperature(:,1),CO2LEM(:,1),'r','LineWidth',1.5);
              subplot(3,2,5);
     p=plot(Temperature(:,1),O2LEM(:,1),'r','LineWidth',1.5);
% plot(Temperature(:,1),species(:,13,1),'b','LineWidth',2);
% plot(Temperature(:,740),species(:,13,740),'r','LineWidth',2);

box on;


%%
% minor Scatter 


hold on;
subplot(3,2,2);
plot(T,CO);
xlabel('Temperature (K)');
 ylabel('Y_{CO}');
subplot(3,2,4);
plot(T,OH);
xlabel('Temperature (K)');
 ylabel('Y_{OH}');
subplot(3,2,6);
plot(T,H2);
xlabel('Temperature (K)');
 ylabel('Y_{H2}');
% subplot(2,1,1);

for i = 1000:40:12500
%      p=plot(position(:,1),Temperature(:,1),'k','LineWidth',2);
     hold on;
%      p=plot(position(:,1),Temperature(:,i),'--r');
% p=plot(species(:,12,1),species(:,11,1),'k','LineWidth',1.5);
% if i<=10000
    subplot(3,2,1);
%      p=plot(Temperature(:,1),species(:,10,1),'r','LineWidth',1.5);
     hold on;
     plot(Temperature(:,i),COLEM(:,i),'k');
     xlabel('Temperature (K)');
 ylabel('Y_{CO}');
         subplot(3,2,3);
%      p=plot(Temperature(:,1),species(:,12,1),'r','LineWidth',1.5);
     hold on;
     plot(Temperature(:,i),OHLEM(:,i),'k');
     xlabel('Temperature (K)');
 ylabel('Y_{OH}');
         subplot(3,2,5);
%      p=plot(Temperature(:,1),species(:,4,1),'r','LineWidth',1.5);
     hold on;
     plot(Temperature(:,i),H2LEM(:,i),'k');
     xlabel('Temperature (K)');
 ylabel('Y_{H2}');
%      p=plot(species(:,12,i),species(:,11,i),'--r');
i
%     
% elseif i>10000
%     subplot(3,2,5);
% %      p=plot(Temperature(:,1),species(:,1,1),'r','LineWidth',1.5);
%      hold on;
%      plot(Temperature(:,i),species2(:,4,i-10000),'k');
% %      xlabel('Temperature (K)');
% %  ylabel('Y_{CO2}');
%          subplot(3,2,3);
% %      p=plot(Temperature(:,1),species(:,11,1),'r','LineWidth',1.5);
%      hold on;
%      plot(Temperature(:,i),species2(:,12,i-10000),'k');
% %      xlabel('Temperature (K)');
% %  ylabel('Y_{CO}');
%          subplot(3,2,1);
% %      p=plot(Temperature(:,1),species(:,5,1),'r','LineWidth',1.5);
%      hold on;
%      plot(Temperature(:,i),species2(:,10,i-10000),'k');
%      xlabel('Y_{CO2}');
%      ylabel('Y_{CO}');
%      title('LEM Results K74');
%      grid on;
     hold on;
%       pause(0.09);
%       delete(p);
% i
% end
 end
    subplot(3,2,1);
     p=plot(Temperature(:,1),COLEM(:,1),'r','LineWidth',1.5);
     box on;
              subplot(3,2,3);
              box on;
     p=plot(Temperature(:,1),OHLEM(:,1),'r','LineWidth',1.5);
              subplot(3,2,5);
              box on;
     p=plot(Temperature(:,1),H2LEM(:,1),'r','LineWidth',1.5);
% plot(Temperature(:,1),species(:,13,1),'b','LineWidth',2);
% plot(Temperature(:,740),species(:,13,740),'r','LineWidth',2);
%%
% subplot(2,1,2);

%      p=plot(position(:,1),Temperature(:,1),'k','LineWidth',1.2);
%  p=plot(Temperature(:,1),species(:,12,1),'k','LineWidth',1.5);
%  hold on;
%      p=plot(position(:,1),Temperature(:,860),'--r');
% y=T(50,10,:);
subplot(3,2,2);
plot(T,CH4,'b');
     xlabel('Temperature (K)');
 ylabel('Y_{CH4}');
hold on;
subplot(3,2,4);
plot(T,CO2,'b');
     xlabel('Temperature (K)');
 ylabel('Y_{CO2}');
hold on;
subplot(3,2,6);
plot(T,O2,'b');
     xlabel('Temperature (K)');
 ylabel('Y_{O2}');
hold on;
%      xlabel('Position (m)');
%      ylabel('Temperature (k)');
     xlabel('Y_{CO2}');
     ylabel('Y_{OH}');
%      title('DNS Data');
%      title('Case K74');


%%
hold on;
dx=0.01/256;
for i=1:256
    
    x(i)=(i-1)*dx;
end
% plot (CO2,CO,'-.b');
T=reshape(T,128,128,256);
y=T(25,68,:);
plot(x,y(:),':r')
ylabel('Temperature (k)');
% plot(x,T);
xlabel('Position (m)');
title('DNS Data');
% grid on;

%%


average=0.5;
variance_normal=0.5;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);



%% 
% plot PDF

cd 'C:\Users\Mojtaba Amini\Desktop';

% load 'PDFs_1.mat';
% load 'PDFs_5.mat';
% load 'PDFs_15.mat';

            for k=1:50
                
                xprime_variable(k)=0.01+(k-1)*(1/50);
                
            end

subplot(3,2,1);
% mean = 0.25 , var = 0.21 
meancount= 13;
varcount=11;
     plot(xprime_variable,PDFs_2(:,meancount,varcount),'k','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_5(:,meancount,varcount),'b','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_15(:,meancount,varcount),'r','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFsDNS(:,meancount,varcount),'-or');
% beta PDF
average=0.25;
variance_normal=0.21;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);
%     plot(xprime_variable,p,'-.b','LineWidth',1.2);

legend ( ' C_{\lambda} = 1.0',  ' C_{\lambda} = 5.0',' C_{\lambda} = 15.0',' DNS');
% legend ( ' LEM',' DNS',' \beta PDF');

subplot(3,2,2);
% mean = 0.25 , var = 0.51 
meancount= 13;
varcount=26;
     plot(xprime_variable,PDFs_2(:,meancount,varcount),'k','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_5(:,meancount,varcount),'b','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_15(:,meancount,varcount),'r','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFsDNS(:,meancount,varcount),'-or');
     average=0.25;
variance_normal=0.51;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);
%     plot(xprime_variable,p,'-.b','LineWidth',1.2);
% legend ( ' C_{\lambda} = 1.0',  ' C_{\lambda} = 5.0',' C_{\lambda} = 15.0',' DNS');
% legend ( ' LEM',' DNS',' Beta');

subplot(3,2,3);
% mean = 0.55 , var = 0.21 
meancount= 23;
varcount=11;
     plot(xprime_variable,PDFs_2(:,meancount,varcount),'k','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_5(:,meancount,varcount),'b','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_15(:,meancount,varcount),'r','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFsDNS(:,meancount,varcount),'-or');    
     average=0.55;
variance_normal=0.21;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);
%     plot(xprime_variable,p,'-.b','LineWidth',1.2);
% legend ( ' C_{\lambda} = 1.0',  ' C_{\lambda} = 5.0',' C_{\lambda} = 15.0',' DNS');
% legend ( ' LEM',' DNS',' Beta');

subplot(3,2,4);
% mean = 0.55 , var = 0.51 
meancount= 23;
varcount=26;
     plot(xprime_variable,PDFs_2(:,meancount,varcount),'k','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_5(:,meancount,varcount),'b','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_15(:,meancount,varcount),'r','LineWidth',1.2);
% hold on;
     plot(xprime_variable,PDFsDNS(:,meancount,varcount),'-or');   
     average=0.55;
variance_normal=0.51;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);
%     plot(xprime_variable,p,'-.b','LineWidth',1.2);
     
% legend ( ' C_{\lambda} = 1.0',  ' C_{\lambda} = 5.0',' C_{\lambda} = 15.0',' DNS');
% legend ( ' LEM',' DNS',' Beta');

subplot(3,2,5);
% mean = 0.75 , var = 0.21 
meancount= 38;
varcount=11;
     plot(xprime_variable,PDFs_2(:,meancount,varcount),'k','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_5(:,meancount,varcount),'b','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_15(:,meancount,varcount),'r','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFsDNS(:,meancount,varcount),'-or');      
     average=0.75;
variance_normal=0.21;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);
%     plot(xprime_variable,p,'-.b','LineWidth',1.2);
% legend ( ' C_{\lambda} = 1.0',  ' C_{\lambda} = 5.0',' C_{\lambda} = 15.0',' DNS');
% legend ( ' LEM',' DNS',' Beta');

subplot(3,2,6);
% mean = 0.75 , var = 0.51 
meancount= 38;
varcount=26;
     plot(xprime_variable,PDFs_2(:,meancount,varcount),'k','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_5(:,meancount,varcount),'b','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFs_15(:,meancount,varcount),'r','LineWidth',1.2);
hold on;
     plot(xprime_variable,PDFsDNS(:,meancount,varcount),'-or');      
     average=0.75;
variance_normal=0.51;
variance=variance_normal*(average)*(1-average);
a=average*(average*(1-average)/variance-1);
b=a/average*(1-average);
p=betapdf(0.01:0.02:0.99,a,b);
%     plot(xprime_variable,p,'-.b','LineWidth',1.2);
% legend ( ' C_{\lambda} = 1.0',  ' C_{\lambda} = 5.0',' C_{\lambda} = 15.0',' DNS');
% legend ( ' LEM',' DNS',' Beta');

subplot(3,2,1);
str='$$\bar{{c^{\prime}}^{2}} = 0.21 $$';
title(str,'Interpreter','latex');
%  ylabel('mean = 0.25');
 str='$$P(c^{*}) \mid \bar{c} = 0.25 $$';
 ylabel(str,'Interpreter','latex');
 subplot(3,2,2);
str='$$\bar{{c^{\prime}}^{2}} = 0.51 $$';
title(str,'Interpreter','latex');
 subplot(3,2,3);
%  ylabel('mean = 0.55');
 str='$$P(c^{*}) \mid \bar{c} = 0.55 $$';
 ylabel(str,'Interpreter','latex');
 subplot(3,2,5);
%  ylabel('mean = 0.55');
 str='$$P(c^{*}) \mid \bar{c} = 0.75 $$';
 ylabel(str,'Interpreter','latex');    
 xlabel('c^{*}');
 subplot(3,2,6);
 xlabel('c^{*}');
 subplot(3,2,2);
 ylim([0 5]);
 subplot(3,2,4);
 ylim([0 5]);
  subplot(3,2,6);
 ylim([0 5]);
%%

maxbin=1600/50;
sum1=zeros(maxbin);
sum12=zeros(maxbin);
sum2=zeros(maxbin);
sum22=zeros(maxbin);
check1=zeros(maxbin,1);
check12=zeros(maxbin,1);
check2=zeros(maxbin,1);
check22=zeros(maxbin,1);

for counter=1:maxbin
    T1=(counter-1)*0.002;
    T2=(counter)*0.002;
    
    for i=1:4000
        for j=1:3850        
            if Temperature(j,i)>=T1 && Temperature(j,i)<=T2

                sum1(counter)=sum1(counter)+species(j,5,i);
                check1(counter)=check1(counter)+1.0;
            end
        end
    end
    bin1(counter)=sum1(counter)/check1(counter);
        for i=1:4000
        for j=1:3850        
            if Temperature(j,i)>=T1 && Temperature(j,i)<=T2

                sum12(counter)=(bin1(counter)-species(j,5,i))^2+sum12(counter);
                check12(counter)=check12(counter)+1.0;
            end
        end
    end
        bin12(counter)=sum12(counter)/(check12(counter)-1);

    for i=4000:6999
        for j=1:3850     
            if Temperature(j,i)>=T1 && Temperature(j,i)<=T2
                sum2(counter)=sum2(counter)+species(j,5,i);
                check2(counter)=check2(counter)+1.0;
            end
        end
     end
%     bin1(counter)=sum/check;
    bin2(counter)=sum2(counter)/check2(counter);
    for i=4000:6999
        for j=1:3850        
            if Temperature(j,i)>=T1 && Temperature(j,i)<=T2
                sum22(counter)=(bin2(counter)-species(j,5,i))^2+sum22(counter);
                check22(counter)=check22(counter)+1.0;
            end
        end
    end
        bin22(counter)=sum22(counter)/(check22(counter)-1);
end


plot(bin1);
hold on;
plot(bin2);
title('Location Average of OH Species');
legend('t=0.04','t=0.07');
figure;
xlabel('Bin');
ylabel('Average');
plot(bin12);
hold on;
plot(bin22);
title('Location Variance of OH Species');
legend('t=0.04','t=0.07');
xlabel('Bin');
ylabel('Variance');

%%


% maxbin=1600/50;
clear -regxp binCO binOH binCO2 binCO2DNS check1 check12 check2 check22;
clear -regxp sum1 sum12 sum2 sum22;
cd 'D:\Selected Ka74 data (1)\Selected Ka74 data';
maxbin=50;
sumOH=zeros(maxbin,1);
checkOH=zeros(maxbin,1);
binOH=zeros(maxbin,1);
sumCO=zeros(maxbin,1);
checkCO=zeros(maxbin,1);
binCO=zeros(maxbin,1);
sumCO2=zeros(maxbin,1);
checkCO2=zeros(maxbin,1);
binCO2=zeros(maxbin,1);

for counter=1:maxbin
    T1=(counter-1)*0.0002;
    T2=(counter)*0.0002;
    
    for i=1:6999
        for j=1:3850        
%             if Temperature(j,i)>=T1 && Temperature(j,i)<=T2
            if position(j,1)>=T1 && position(j,1)<=T2
                sumOH(counter,1)=sumOH(counter,1)+species(j,5,i);
                sumCO(counter,1)=sumCO(counter,1)+species(j,11,i);
                sumCO2(counter,1)=sumCO2(counter,1)+species(j,12,i);
                
                checkOH(counter,1)=checkOH(counter,1)+1.0;
                checkCO(counter,1)=checkCO(counter,1)+1.0;
                checkCO2(counter,1)=checkCO2(counter,1)+1.0;
            end
        end
    end
    binOH(counter,1)=sumOH(counter,1)/checkOH(counter,1);
    binCO(counter,1)=sumCO(counter,1)/checkCO(counter,1);
    binCO2(counter,1)=sumCO2(counter,1)/checkCO2(counter,1);

end

CO2=reshape(CO2,128,128,256);
CO=reshape(CO,128,128,256);
OH=reshape(OH,128,128,256);
% y=T(25,68,:);
% plot(x,y(:),':r')


maxbin=50;
sumOHDNS=zeros(maxbin,1);
checkOHDNS=zeros(maxbin,1);
binOHDNS=zeros(maxbin,1);
sumCODNS=zeros(maxbin,1);
checkCODNS=zeros(maxbin,1);
binCODNS=zeros(maxbin,1);
sumCO2DNS=zeros(maxbin,1);
checkCO2DNS=zeros(maxbin,1);
binCO2DNS=zeros(maxbin,1);


for counter=1:50
    T1=(counter-1)*5;
    T2=(counter)*5;
    
    for i=1:128
        for j=1:128
            for k=3:253
                
                if k>=T1 && k<=T2
                sumOHDNS(counter,1)=sumOHDNS(counter,1)+OH(i,j,k);
                sumCODNS(counter)=sumCODNS(counter,1)+CO(i,j,k);
                sumCO2DNS(counter,1)=sumCO2DNS(counter,1)+CO2(i,j,k);
                
                checkOHDNS(counter,1)=checkOHDNS(counter,1)+1.0;
                checkCODNS(counter,1)=checkCODNS(counter,1)+1.0;
                checkCO2DNS(counter,1)=checkCO2DNS(counter,1)+1.0;
                    
                end
                
            end
        end
    end
    binOHDNS(counter,1)=sumOHDNS(counter,1)/checkOHDNS(counter,1);
    binCODNS(counter,1)=sumCODNS(counter,1)/checkCODNS(counter,1);
    binCO2DNS(counter,1)=sumCO2DNS(counter,1)/checkCO2DNS(counter,1);
end
              
        






plot(binOH);
hold on;
plot(binOHDNS);
title('Spatial Average of OH Species');
xlabel('Bin');
ylabel('Average');
legend('LEM','DNS');
figure;
plot(binCO);
hold on;
plot(binCODNS);
title('Spatial Average of CO Species');
xlabel('Bin');
ylabel('Average');
legend('LEM','DNS');
figure;
plot(binCO2);
hold on;
plot(binCO2DNS);
title('Spatial Average of CO2 Species');
xlabel('Bin');
ylabel('Average');
legend('LEM','DNS');

%%
% rea .txt file for drawing f(l)
clear maximum;
clear minimum;
clear lengthPDF;
clear l;
cd 'F:\test\DNS_cases\LaminarFlameSpeedNewMechanism\LEMNewFrequencyEquation5';
fileID = fopen('Size_TripleMap.txt','r');
formatSpec = '%f';
tripletSizes = fscanf(fileID,formatSpec);
lengharray=length(tripletSizes);
maximum=max(tripletSizes);
minimum=min(tripletSizes);
lengthPDF=zeros((301),1);
k=1;
for i=48:348
    for j=1:lengharray
        if tripletSizes(j)==i
            lengthPDF(k)=lengthPDF(k)+1; 
        end
    end
    k=k+1;
end
%%
xlint= 0.1196; %cm
Ret=32;
Neta=1.08;
xlk=Neta* xlint/(Ret^(3/4));
%%
l=linspace(xlk,xlint,301);
bar(l,lengthPDF);
% title('f(l) Data');
% Dx=0.01/3850;
% PDFA=((xlint)^(5/3)*(xlk)^(-5/3))/((xlint/xlk)^(5/3)-1.0);
% PDFB=-((xlint)^(5/3))/((xlint/xlk)^(5/3)-1.0);
% zp=linspace(0,1,100);
% for i=1:100
% min(i)=(((zp(i)-PDFA)/PDFB)^(-3/5))/(Dx*100);
% end
% max=((1-PDFA/PDFB)^(-3/5))/(Dx*100);

% l=linspace(xlk,xlint,maximum-minimum+1);
for i=1:301
PDF(1,i)=((5/3)*l(i)^(-8/3))/(xlk^(-5/3)-xlint^(-5/3));

end
hold on;
plot(l,PDF,'r');
title('f(l) Equation VS Numerical Data');
xlabel('Eddy Size');
ylabel('Number of Tripletmap');

ylim([0 340]);

%%

% PDF
c=0.01:0.02:0.99;
for i=1:8
    for j=1:8
        for k=1:16
            average(i,j,k)= sum(c'.*p(:,i,j,k)/50);
            variance(i,j,k)= sum(p(:,i,j,k)/50.*(c'-average(i,j,k)).^2);
        end
    end
end
%%
p=50*samplesize/sum(samplesize(:,1,1,1));
plot(0:1/49:1,p(:,5,6,6))



%%

phi=[0.4,0.6,0.8,1.0,1.2,1.4,1.6];
GRI=[0.008351,0.121569,0.299739,0.406322,0.352709,0.145654,0.082012];
 plot(phi,GRI,'--k*');

hold on;
xlabel('Equivalence Ratio \phi');
ylabel('Laminar Flame Speed S_{L} (m/s)');

SmookeCantera=[0.013236,0.160884,0.394272,0.611008,0.703705,0.573016,0.314136];
plot(phi,SmookeCantera,'--rs');

hold on;
box on;
smookeME=[0.01270656,0.158842,0.39032928,0.598859,0.6825938,0.561555,0.3047119];
title('Laminar Flame Speed');
% plot(phi,smookeME,'--bo');

legend('GRI','SMOOKE');


%%

phi=[0.4,0.6,0.8,1.0,1.2,1.4,1.6];
GRI=[0.008351,0.121569,0.299739,0.406322,0.352709,0.145654,0.082012];
plot(phi,GRI,'--rs');

hold on;
xlabel('Equivalence Ratio \phi');
ylabel('Laminar Flame Speed S_{L} (m/s)');

hold on;

results=[0.00793354,0.11913762,0.284752,0.394132,0.33154646,0.142740,0.078731];
title('Laminar Flame Speed');
plot(phi,results,'--bo');

legend('GRI','Results');


%%

L=0.092;
Re=120;
Neta= 20.0;
kolmogrove=Neta*L*Re^(-3/4);
Clanda=15.0;
kineticViscosity=0.1818;
freq=(54/5)*((kineticViscosity*Re)/(Clanda*L^3))*(((L/kolmogrove)^(5/3)-1)/(1-(kolmogrove/L)^(4/3)));
%%
figure;

first=[1,2,5,10,15,20];
firsty=[1.2088e+8,3.8355e+7,8.5246e+6,2.7868e+06,1.4668e+06,9.3670e+05];
semilogy(first,firsty,'.-r');
hold on;
second=[1,2,5,10,15,20];
secondy=[2.4176e+07,7.6711e+6,1.7049e+6,5.5737e+5,2.9337e+05,1.8734e+05];
semilogy(second,secondy,'--k');
hold on;
third=[1,2,5,10,15,20];
thirdy=[8.0586e+06,2.5570e+6,5.6831e+5,1.8579e+5,9.7789e+4,6.2446e+4];
semilogy(third,thirdy,'--bs');
grid on;
legend('C_{\lambda} = 1.0','C_{\lambda} = 5.0','C_{\lambda} = 15');
ax1 = gca; % current axes
% ax1.XColor = 'r';
xlabel('N_{\eta}');
ylabel ('\lambda Frequency');

%%
% create file for PDF

% cd 'F:\test\DNS_cases\LaminarFlameSpeedNewMechanism\LEMNewFrequencyEquation5\Turbulence_K74';

% Number of Profile saved


% read all data of ember code from h5py file


% write Progress variable in a csv file
% using LEM data to create PDF
% fid = fopen('F:\test\DNS_cases\LaminarFlameSpeedNewMechanism\LEMNewFrequencyEquation8 (copy)\interpedc.csv','w');
% for j=1:6999
%     for i=1:3850
%         fprintf(fid,'%f, ',species(i,13,j));
%         fprintf(fid,'\n'); 
%     end
%    
% end
% fclose(fid);
fid = fopen('F:\test\DNS_cases\LaminarFlameSpeedNewMechanism\LEMNewFrequencyEquation8 (copy)\x.csv','w');

    for i=1:3850
        fprintf(fido,'%f, ',position(i,1));
        
        fprintf(fido,'\n'); 
    end
fclose(fido);

%%

% Post processing of PDF Data

cd 'D:\Create PDF & SDR';

% load PDF.csv;
data=csvread('Only_PDF.csv');
C_tild=csvread('Only_Ctild.csv');
C_var=csvread('Only_Cvar.csv');



hold on;
mean=[0.11 0.31 0.51 0.71 0.91];
var=[0.25 0.25 0.25 0.25 0.25];


for count=1:5
    
    for i=1:53:132500

        if mean(count)==data(i,1) && var(count)==data(i+1,1)
            
            
            for k=1:50
                
                xprime_variable(k)=0.01+(k-1)*(1/50);
                
            end
            x_variable(1)=(1/50)/2;
            x_variable(51)=49*(1/50)+(1/50)/2;
            x_variable=0:0.02:1;
            for j=1:51
                desire_PDF(j)=50*data(i+1+j,1);
%               sum(desire_PDF);
            end
            plot_PDF(1)=desire_PDF(1)+desire_PDF(1+1)/2;
            plot_PDF(50)=desire_PDF(50)+desire_PDF(50+1)/2;
            for j=2:49
                plot_PDF(j)=(desire_PDF(50)+desire_PDF(j+1))/2;
%               sum(desire_PDF);
            end
            break;
        end
    end
    count
%     sum(desire_PDF)
   subplot(5,3,3*count-2);
      hold on;
     plot(xprime_variable,LEM_PDFs(:,10*(count-1)+6,13),'-.k');
%       a=mean(count)*(mean(count)*(1-mean(count))/var(count)-1);
%       b=a/mean(count)*(1-mean(count));
%       p=betapdf(xprime_variable,a,b);
%       hold on;
%       plot(xprime_variable,p,'r');

                  plot(xprime_variable,PDFs(:,10*(count-1)+6,13),'k');
%                   legend('LEM PDF','DNS PDF');
%       pause(0.5);
%       delete(p1);
%       delete(p2);
    %     xlim([0 1]);
%     ylim([0 5]);
    str='$$P(c^{*}) \mid \bar{c}= $$';
    ylabel(str,'Interpreter','latex');
%     if count==1
    str='$$\bar{{c^{\prime}}^{2}} $$';
    title(str,'Interpreter','latex');
%     end
%     if count==5
%         
        xlabel('c^{*}');
%     end
%     ylabel(num2str(mean(count)));
    
end




%% 2 for test
% Post processing of PDF Data

cd 'D:\Create PDF & SDR';

% load PDF.csv;
data=csvread('Only_PDF.csv');
C_tild=csvread('Only_Ctild.csv');
C_var=csvread('Only_Cvar.csv');

hold on;
mean=0.87;
var=0.27;

% for count=1:50
    
    for i=1:53:132500

        if mean==data(i,1) && var==data(i+1,1)
            
            
            for k=1:49
                
                xprime_variable(k+1)=k*(1/50);
                
            end
            x_variable(1)=(1/50)/2;
            x_variable(51)=49*(1/50)+(1/50)/2;
            x_variable=0:0.02:1;
            for j=1:51
                desire_PDF(j)=50*data(i+1+j,1);
%               sum(desire_PDF);
            end
            break;
        end
    end
%     count
%     sum(desire_PDF)
   subplot(3,3,7);
      hold on;
      p1=plot(x_variable,desire_PDF,'r');
      
      p2=plot(xprime_variable,PDFs(:,43,13),'b--');
%       pause(1.5);
%       delete(p1);
%       delete(p2);
    %     xlim([0 1]);
    ylim([0 20]);
%     str='$$P(c^{*}) \mid \bar{c}= $$';
%     ylabel(str,'Interpreter','latex');
%     if count==1
%     str='$$\bar{{c^{\prime}}^{2}} $$';
%     title(str,'Interpreter','latex');
%     end
%     if count==5
%         
%         xlabel('c^{*}');
%     end
%     ylabel(num2str(mean(count)));
    
% end

%%


count=0.01:0.01:1;

for k=1:100


    for i=1:53:132500

        if 0.2==data(i,1) && var(count(k))==data(i+1,1)
                    
            for k=1:49    
                xprime_variable(k+1)=k*(1/50);
            end
            x_variable(1)=(1/50)/2;
            x_variable(51)=49*(1/50)+(1/50)/2;
            x_variable=0:0.02:1;
            for j=1:51
                desire_PDF(j,k)=50*data(i+1+j,1);
%               sum(desire_PDF);
            end
            break;
        end
    end
end
    for k=1:100
    
   p= plot(x_variable,desire_PDF,'--r');
    
    hold on;
      pause(0.1);
      delete(p);
    



end



%%

gradc=csvread('gradc.csv');
for i=1:(length(gradc)/2)
    
   x(i)=gradc(2*i-1,1); 
   y(i)=gradc(2*i,1); 
end

plot(x,y);




%%

x= 1:31;
y = 3*x;
subplot(3,2,[1,3,5]);
plot(x,y,'-rs');
hold on;
x= 8:23;
y = 3*x;
plot(x,y,'-ks');
hold on;
x=8;
y=24;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=13;
y=39;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=18;
y=54;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=23;
y=69;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
%---------

set(gca,'XTick',[],'YTick',[]); % this works okay
ylim([-10,100])
xlim([-1,34]);
str = {'Initial Domain'};
text(4,90,str)

%--------
subplot(3,2,2);
x= 1:31;
y = 3*x;
plot(x,y,'-rs');
hold on;
x= 13:23;
y = 3*x;
plot(x,y,'-ks');
hold on;
x=8;
y=24;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=13;
y=39;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=18;
y=54;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=23;
y=69;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 


set(gca,'XTick',[],'YTick',[]); % this works okay
ylim([-10,100])
xlim([-1,34]);
hold on;
x=8:13;
y=9*x-6*8;
plot(x,y,'-bo');
str = {'\fontsize{8} X_{0}\leq X \leq X_{0}+ l/3'};
text(17,10,str)
str = {'\fontsize{16} 1'};
text(1,75,str)

%-----------

subplot(3,2,4);
x= 1:31;
y = 3*x;
plot(x,y,'-rs');
hold on;
x= 18:23;
y = 3*x;
plot(x,y,'-ks');
hold on;
x=8;
y=24;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=13;
y=39;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=18;
y=54;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=23;
y=69;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 


set(gca,'XTick',[],'YTick',[]); % this works okay
ylim([-10,100])
xlim([-1,34]);
hold on;
x=8:13;
y=9*x-6*8;
plot(x,y,'-bo');
hold on;
x=13:18;
y=-9*x+12*8+6*15;
plot(x,y,'-bo');
str = {'\fontsize{8} X_{0}+ l/3\leq X \leq X_{0}+ 2l/3'};
text(13,6,str)
str = {'\fontsize{16} 2'};
text(1,75,str)

%--------

subplot(3,2,6);
x= 1:31;
y = 3*x;
plot(x,y,'-rs');
hold on;
% x= 18:23;
% y = 3*x;
% plot(x,y,'-ks');
% hold on;
x=8;
y=24;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=13;
y=39;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=18;
y=54;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 
hold on;
x=23;
y=69;
z=plot(x,y,'-bo');
set(z,'MarkerFaceColor', get(z,'Color')); 


set(gca,'XTick',[],'YTick',[]); % this works okay
ylim([-10,100])
xlim([-1,34]);
hold on;
x=8:13;
y=9*x-6*8;
plot(x,y,'-bo');
hold on;
x=13:18;
y=-9*x+12*8+6*15;
plot(x,y,'-bo');
hold on;
x=18:23;
y=9*x-6*8-6*15;
plot(x,y,'-bo');

str = {'\fontsize{8} X_{0}+ 2l/3\leq X \leq X_{0}+ l'};
text(13,6,str)
str = {'\fontsize{16} 3'};
text(1,75,str)

%%


x=[0.1 10];
y=[10 0.1];
loglog(x, y,'-.k','linewidth',1.0)
hold on;
x=[1 2*10e+2];
y=[1 1];
loglog(x, y,'-.k','linewidth',1.0)
hold on;
x=[1 2*10e+2];
y=[1 2*10e+2];
loglog(x, y,'-.k','linewidth',1.0)
hold on;
x=[0.1 2*10e+2];
y=[10 2*10e+1];
loglog(x, y,'-.k','linewidth',1.0)
hold on;
x=[1 2*10e+2];
y=[1 10];
loglog(x, y,'-.k','linewidth',1.0)
hold on;
clear x;
clear y;
x=1.3;
y=3.7;
p = plot(x,y,'rs');
set(p,'MarkerFaceColor', get(p,'Color'));
str = {'A1'};
text(1.5,3.8,str)
str = {'A2'};
text(1.8,3.8,str)
str = {'B'};
text(2.1,3.8,str)
str = {'C'};
text(3,3.8,str)
str = {'D'};
text(1.2,19,str)
x=1;
y=18;
p = plot(x,y,'rs');
set(p,'MarkerFaceColor', get(p,'Color'));

x=0.48;
y=210;
p = plot(x,y,'bs');
set(p,'MarkerFaceColor', get(p,'Color'));

x=2.4;
y=24;
p = plot(x,y,'bs');
set(p,'MarkerFaceColor', get(p,'Color'));

x=1.2;
y=13;
p = plot(x,y,'bs');
set(p,'MarkerFaceColor', get(p,'Color'));
str = {'A3'};
text(0.55,212,str)
str = {'Thickened Flame'};
text(0.38,812,str)
str = {'Thickened-Wrinkled Flame'};
% dim = [.74 .56 .1 .1];
% annotation('textbox',dim,'String',str,'FaceColor','blue','FitBoxToText','on');
text(12,18,str)
str = {'Wrinkled Flamelet'};
text(60,1,str)
str = {'Laminar Combustion'};
text(0.13,0.4,str)
x=1.0;
y=66;
p = plot(x,y,'bs');
set(p,'MarkerFaceColor', get(p,'Color'));
% set(gca,'Color',[200 200 200])
ylim([10e-2,2*10e+2]);
xlim([10e-2,2*10e+2]);
hold on;
xlabel('L / \delta_{f}');
ylabel('u^{\prime} / s_{L}');




%%

clc;
xlk=0.0097;
nu=0.1562;
re=32;
ep=(nu^3)/(xlk^4);
ta=sqrt(nu/ep);

ta*re



%%

% Scalar Dissipation rate calculation
clc;
% based on CO2 mass fraction
cd 'F:\Thesis\DNS data\K75_k1';

load 'Yi[14].mat';
%%
CO2=reshape(CO2,128,128,256);


%%

clc;
twodx=position(4,1)-position(2,1);
dif=zeros(1031,12500);
for j=4000:12500
    
   for i=2:1031
       
       dif(i,j)=((CO2LEM(i+1,j)-CO2LEM(i-1,j))/twodx)^2;
                
   end
       
end

twodxDNS=0.01/255;
difDNS=zeros(126,126,254);

for i=2:254
    
   for j=2:126
       
       for k=2:126
           
           xdiff=(CO2(k,j,i+1)-CO2(k,j,i-1))/twodxDNS;
           ydiff=(CO2(k,j+1,i)-CO2(k,j-1,i))/twodxDNS;
           zdiff=(CO2(k+1,j,i)-CO2(k-1,j,i))/twodxDNS;
           difDNS(k,j,i)=xdiff^2+ydiff^2+zdiff^2;
           
       end
                
   end
       
end

for i=4000:12500
    
%     if i<=10000
    maxInCO2total(i)=max(CO2LEM(:,i));
    minInCO2total(i)=min(CO2LEM(:,i));
%     elseif i>10000
%     maxInCO2total(i)=max(species2(:,12,i-10000));
%     minInCO2total(i)=min(species2(:,12,i-10000));   
%     end
end

maxCO2Bin=max(maxInCO2total(:));
minCO2Bin=min(minInCO2total(:));

numbin=100;
binsize=(maxCO2Bin-minCO2Bin)/numbin;
sumSDR=zeros(numbin,1);
checkSDR=zeros(numbin,1);
binSDR=zeros(numbin,1);

sumSDRDNS=zeros(numbin,1);
checkSDRDNS=zeros(numbin,1);
binSDRDNS=zeros(numbin,1);

for counter=1:numbin
    F1=(counter-1)*binsize;
    F2=(counter)*binsize;
    
    for i=4000:12500
        for j=1:1031
            if CO2LEM(j,i)>=F1 && CO2LEM(j,i)<=F2

                sumSDR(counter,1)=sumSDR(counter,1)+dif(j,i);
                
                xseiesdata(counter)=((F1+F2)/2)/maxCO2Bin;

                checkSDR(counter,1)=checkSDR(counter,1)+1.0;
            end
        end
    end
    
    binSDR(counter,1)=sumSDR(counter,1)/checkSDR(counter,1);
end


for counter=1:numbin
    F1=(counter-1)*binsize;
    F2=(counter)*binsize;
    
    for k=1:126
        for j=1:126
            for i=1:254
                
                if CO2(k,j,i)>=F1 && CO2(k,j,i)<=F2
                sumSDRDNS(counter,1)=sumSDRDNS(counter,1)+difDNS(k,j,i);


                checkSDRDNS(counter,1)=checkSDRDNS(counter,1)+1.0;
                end
                
            end
        end
    end
    binSDRDNS(counter,1)=sumSDRDNS(counter,1)/checkSDRDNS(counter,1);
end

% maxLEM=max(binSDR);
% maxDNS=max(binSDRDNS);
% binSDR=binSDR/maxLEM;
% binSDRDNS=binSDRDNS/maxDNS;

plot(xseiesdata,binSDR,'r');
hold on;
plot(xseiesdata,binSDRDNS,'b');

ylabel('\chi_{c}|c^{*}');
xlabel('c^{*}');
title('Case A2');
legend('DNS SDR','LEM SDR');

%%
plot(sumwdotCO2);
xlabel('Sample Number')
    str='$$\int\rho\dot{\omega}_{CO2}$$';
    ylabel(str,'Interpreter','latex','FontSize',15);




