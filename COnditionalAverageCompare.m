
clc;

for i=1:19999
    
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

clear maxInCO2total;
clear minInCO2total;
clear i;

%%

% c lambda = 1

numbin=100;
binsize=(maxCO2Bin-minCO2Bin)/numbin;
sumOH=zeros(numbin,1);
checkOH=zeros(numbin,1);
binOH=zeros(numbin,1);
sumCO=zeros(numbin,1);
checkCO=zeros(numbin,1);
binCO=zeros(numbin,1);
sumCO2=zeros(numbin,1);
checkCO2=zeros(numbin,1);
binCO2=zeros(numbin,1);
sumO2=zeros(numbin,1);
checkO2=zeros(numbin,1);
binO2=zeros(numbin,1);
sumH2=zeros(numbin,1);
checkH2=zeros(numbin,1);
binCH4=zeros(numbin,1);
sumCH4=zeros(numbin,1);
checkCH4=zeros(numbin,1);
binH2=zeros(numbin,1);
sumTemp=zeros(numbin,1);
checkTemp=zeros(numbin,1);
binTemp=zeros(numbin,1);
xseiesdata=zeros(numbin,1);

for counter=1:numbin
    F1=(counter-1)*binsize;
    F2=(counter)*binsize;
    
    for i=1000:12500
        for j=1:1032   
            if CO2LEM(j,i)>=F1 && CO2LEM(j,i)<=F2
            
                sumOH(counter,1)=sumOH(counter,1)+OHLEM(j,i);
                  sumCO(counter,1)=sumCO(counter,1)+COLEM(j,i);
                  sumCO2(counter,1)=sumCO2(counter,1)+CO2LEM(j,i);
                  sumO2(counter,1)=sumO2(counter,1)+O2LEM(j,i);
                sumH2(counter,1)=sumH2(counter,1)+H2LEM(j,i);
                  sumCH4(counter,1)=sumCH4(counter,1)+CH4LEM(j,i);
                sumTemp(counter,1)=sumTemp(counter,1)+Temperature(j,i);

                xseiesdata(counter)=((F1+F2)/2)/maxCO2Bin;

                checkOH(counter,1)=checkOH(counter,1)+1.0;
                checkCO(counter,1)=checkCO(counter,1)+1.0;
                checkCO2(counter,1)=checkCO2(counter,1)+1.0;
                  checkO2(counter,1)=checkO2(counter,1)+1.0;
                  checkH2(counter,1)=checkH2(counter,1)+1.0;
                  checkCH4(counter,1)=checkCH4(counter,1)+1.0;
                  checkTemp(counter,1)=checkTemp(counter,1)+1.0;

            end
        end
    end
    binOH(counter,1)=sumOH(counter,1)/checkOH(counter,1);
    binCO(counter,1)=sumCO(counter,1)/checkCO(counter,1);
    binCO2(counter,1)=sumCO2(counter,1)/checkCO2(counter,1);
    binO2(counter,1)=sumO2(counter,1)/checkO2(counter,1);
    binH2(counter,1)=sumH2(counter,1)/checkH2(counter,1);
    binCH4(counter,1)=sumCH4(counter,1)/checkCH4(counter,1);
    binTemp(counter,1)=sumTemp(counter,1)/checkTemp(counter,1);


end


%%

% c lambda = 5
% new data

numbin=100;
binsize=(maxCO2Bin-minCO2Bin)/numbin;
sumOHnew=zeros(numbin,1);
checkOHnew=zeros(numbin,1);
binOHnew=zeros(numbin,1);
sumCOnew=zeros(numbin,1);
checkCOnew=zeros(numbin,1);
binCOnew=zeros(numbin,1);
sumCO2new=zeros(numbin,1);
checkCO2new=zeros(numbin,1);
binCO2new=zeros(numbin,1);
sumO2new=zeros(numbin,1);
checkO2new=zeros(numbin,1);
binO2new=zeros(numbin,1);
sumH2new=zeros(numbin,1);
checkH2new=zeros(numbin,1);
binCH4new=zeros(numbin,1);
sumCH4new=zeros(numbin,1);
checkCH4new=zeros(numbin,1);
binH2new=zeros(numbin,1);
sumTempnew=zeros(numbin,1);
checkTempnew=zeros(numbin,1);
binTempnew=zeros(numbin,1);
xseiesdata=zeros(numbin,1);

for counter=1:numbin
    F1=(counter-1)*binsize;
    F2=(counter)*binsize;
    
    for i=1:12139
        for j=1:3850        
            if species(j,12,i)>=F1 && species(j,12,i)<=F2
            
                sumOHnew(counter,1)=sumOHnew(counter,1)+species(j,5,i);
                  sumCOnew(counter,1)=sumCOnew(counter,1)+species(j,11,i);
                  sumCO2new(counter,1)=sumCO2new(counter,1)+species(j,12,i);
                  sumO2new(counter,1)=sumO2new(counter,1)+species(j,4,i);
                sumH2new(counter,1)=sumH2new(counter,1)+species(j,1,i);
                  sumCH4new(counter,1)=sumCH4new(counter,1)+species(j,10,i);
                sumTempnew(counter,1)=sumTempnew(counter,1)+Temperature(j,i);

                xseiesdata(counter)=((F1+F2)/2)/maxCO2Bin;

                checkOHnew(counter,1)=checkOHnew(counter,1)+1.0;
                checkCOnew(counter,1)=checkCOnew(counter,1)+1.0;
                checkCO2new(counter,1)=checkCO2new(counter,1)+1.0;
                  checkO2new(counter,1)=checkO2new(counter,1)+1.0;
                  checkH2new(counter,1)=checkH2new(counter,1)+1.0;
                  checkCH4new(counter,1)=checkCH4new(counter,1)+1.0;
                  checkTempnew(counter,1)=checkTempnew(counter,1)+1.0;

            end
        end
    end
    binOHnew(counter,1)=sumOHnew(counter,1)/checkOHnew(counter,1);
    binCOnew(counter,1)=sumCOnew(counter,1)/checkCOnew(counter,1);
    binCO2new(counter,1)=sumCO2new(counter,1)/checkCO2new(counter,1);
    binO2new(counter,1)=sumO2new(counter,1)/checkO2new(counter,1);
    binH2new(counter,1)=sumH2new(counter,1)/checkH2new(counter,1);
    binCH4new(counter,1)=sumCH4new(counter,1)/checkCH4new(counter,1);
    binTempnew(counter,1)=sumTempnew(counter,1)/checkTempnew(counter,1);


end

%%


% c lambda =15
% new data

numbin=100;
binsize=(maxCO2Bin-minCO2Bin)/numbin;
sumOHnew1=zeros(numbin,1);
checkOHnew1=zeros(numbin,1);
binOHnew1=zeros(numbin,1);
sumCOnew1=zeros(numbin,1);
checkCOnew1=zeros(numbin,1);
binCOnew1=zeros(numbin,1);
sumCO2new1=zeros(numbin,1);
checkCO2new1=zeros(numbin,1);
binCO2new1=zeros(numbin,1);
sumO2new1=zeros(numbin,1);
checkO2new1=zeros(numbin,1);
binO2new1=zeros(numbin,1);
sumH2new1=zeros(numbin,1);
checkH2new1=zeros(numbin,1);
binCH4new1=zeros(numbin,1);
sumCH4new1=zeros(numbin,1);
checkCH4new1=zeros(numbin,1);
binH2new1=zeros(numbin,1);
sumTempnew1=zeros(numbin,1);
checkTempnew1=zeros(numbin,1);
binTempnew1=zeros(numbin,1);
xseiesdata=zeros(numbin,1);

for counter=1:numbin
    F1=(counter-1)*binsize;
    F2=(counter)*binsize;
    
    for i=1:20000
        for j=1:3850        
            if i<=10000
                if species(j,12,i)>=F1 && species(j,12,i)<=F2
            
                  sumOHnew1(counter,1)=sumOHnew1(counter,1)+species(j,5,i);
                  sumCOnew1(counter,1)=sumCOnew1(counter,1)+species(j,11,i);
                  sumCO2new1(counter,1)=sumCO2new1(counter,1)+species(j,12,i);
                  sumO2new1(counter,1)=sumO2new1(counter,1)+species(j,4,i);
                  sumH2new1(counter,1)=sumH2new1(counter,1)+species(j,1,i);
                  sumCH4new1(counter,1)=sumCH4new1(counter,1)+species(j,10,i);
                  sumTempnew1(counter,1)=sumTempnew1(counter,1)+Temperature(j,i);               

                xseiesdata(counter)=((F1+F2)/2)/maxCO2Bin;

                checkOHnew1(counter,1)=checkOHnew1(counter,1)+1.0;
                checkCOnew1(counter,1)=checkCOnew1(counter,1)+1.0;
                checkCO2new1(counter,1)=checkCO2new1(counter,1)+1.0;
                  checkO2new1(counter,1)=checkO2new1(counter,1)+1.0;
                  checkH2new1(counter,1)=checkH2new1(counter,1)+1.0;
                  checkCH4new1(counter,1)=checkCH4new1(counter,1)+1.0;
                  checkTempnew1(counter,1)=checkTempnew1(counter,1)+1.0;
                end
            end
            if i>10000
                if species2(j,12,i-10000)>=F1 && species2(j,12,i-10000)<=F2
            
                  sumOHnew1(counter,1)=sumOHnew1(counter,1)+species2(j,5,i-10000);
                  sumCOnew1(counter,1)=sumCOnew1(counter,1)+species2(j,11,i-10000);
                  sumCO2new1(counter,1)=sumCO2new1(counter,1)+species2(j,12,i-10000);
                  sumO2new1(counter,1)=sumO2new1(counter,1)+species2(j,4,i-10000);
                  sumH2new1(counter,1)=sumH2new1(counter,1)+species2(j,1,i-10000);
                  sumCH4new1(counter,1)=sumCH4new1(counter,1)+species2(j,10,i-10000);
                  sumTempnew1(counter,1)=sumTempnew1(counter,1)+Temperature(j,i-10000);               

                xseiesdata(counter)=((F1+F2)/2)/maxCO2Bin;

                checkOHnew1(counter,1)=checkOHnew1(counter,1)+1.0;
                checkCOnew1(counter,1)=checkCOnew1(counter,1)+1.0;
                checkCO2new1(counter,1)=checkCO2new1(counter,1)+1.0;
                  checkO2new1(counter,1)=checkO2new1(counter,1)+1.0;
                  checkH2new1(counter,1)=checkH2new1(counter,1)+1.0;
                  checkCH4new1(counter,1)=checkCH4new1(counter,1)+1.0;
                  checkTempnew1(counter,1)=checkTempnew1(counter,1)+1.0;
                end
            end
        end
    end
    binOHnew1(counter,1)=sumOHnew1(counter,1)/checkOHnew1(counter,1);
    binCOnew1(counter,1)=sumCOnew1(counter,1)/checkCOnew1(counter,1);
    binCO2new1(counter,1)=sumCO2new1(counter,1)/checkCO2new1(counter,1);
    binO2new1(counter,1)=sumO2new1(counter,1)/checkO2new1(counter,1);
    binH2new1(counter,1)=sumH2new1(counter,1)/checkH2new1(counter,1);
    binCH4new1(counter,1)=sumCH4new1(counter,1)/checkCH4new1(counter,1);
    binTempnew1(counter,1)=sumTempnew1(counter,1)/checkTempnew1(counter,1);


end

%%
figure;
subplot(3,1,2);
plot(xseiesdata,binOH,'b');
hold on;
% plot(xseiesdata,binOHnew,'-.k');
% hold on;
% plot(xseiesdata,binOHnew1,'b');
% hold on;
plot(xseiesdata,binOHDNS,'r');
title('Conditional Average of OH Species');
xlabel('Normalized CO2 Mass Fraction');
ylabel('OH Mass Fraction Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
legend('LEM','DNS');

subplot(3,1,1);
% figure;
plot(xseiesdata,binCO,'k');
hold on;
% plot(xseiesdata,binCOnew,'-.k');
% hold on;
% plot(xseiesdata,binCOnew1,'b');
% hold on;
plot(xseiesdata,binCODNS,'r');
title('Conditional Average of CO Species');
xlabel('Normalized CO2 Mass Fraction');
ylabel('CO Mass Fraction Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
legend('LEM','DNS');
subplot(3,1,3);
% figure;
plot(xseiesdata,binH2,'k');
hold on;
% plot(xseiesdata,binH2new,'-.k');
% hold on;
% plot(xseiesdata,binH2new1,'b');
% hold on;
plot(xseiesdata,binH2DNS,'r');
title('Conditional Average of H2 Species');
xlabel('Normalized CO2 Mass Fraction');
ylabel('H2 Mass Fraction Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
legend('LEM','DNS');

% figure;
% plot(xseiesdata,binCO2,'k');
% hold on;
% plot(xseiesdata,binCO2new,'-.k');
% hold on;
% plot(xseiesdata,binCO2new1,'b');
% hold on;
% plot(xseiesdata,binCO2DNS,'r');
% title('Conditional Average of CO2 Species');
% xlabel('Normalized CO2 Mass Fraction');
% ylabel('CO2 Mass Fraction Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
figure;
subplot(3,1,3);
plot(xseiesdata,binO2,'k');
hold on;
% plot(xseiesdata,binO2new,'-.k');
% hold on;
% plot(xseiesdata,binO2new1,'b');
% hold on;
plot(xseiesdata,binO2DNS,'r');
title('Conditional Average of O2 Species');
xlabel('Normalized CO2 Mass Fraction');
ylabel('O2 Mass Fraction Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
legend('LEM','DNS');
subplot(3,1,2);
% figure;
plot(xseiesdata,binCH4,'k');
hold on;
% plot(xseiesdata,binCH4new,'-.k');
% hold on;
% plot(xseiesdata,binCH4new1,'b');
% hold on;
plot(xseiesdata,binCH4DNS,'r');
title('Conditional Average of CH4 Species');
xlabel('Normalized CO2 Mass Fraction');
ylabel('CH4 Mass Fraction Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
legend('LEM','DNS');
subplot(3,1,1);
% figure;
plot(xseiesdata,binTemp,'k');
hold on;
% plot(xseiesdata,binTempnew,'-.k');
% hold on;
% plot(xseiesdata,binTempnew1,'b');
% hold on;
plot(xseiesdata,binTempDNS,'r');
title('Conditional Average of Temperature ');
xlabel('Normalized CO2 Mass Fraction');
ylabel('Temperature Average');
% legend('LEM C_{\lambda} = 1.0','LEM C_{\lambda} = 5.0','LEM C_{\lambda} = 15.0','DNS');
legend('LEM','DNS');

%%
clear -regxp binOH sumOH checkOH binCO sumCO checkCO binCO2 sumCO2 checkCO2;
clear -regxp binO2 sumO2 checkO2 binH2 sumH2 checkH2 binCH4 sumCH4 checkCH4;
clear -regxp binTemp sumTemp checkTemp binsize numbin counter F1 F2 i j;

%%

cd 'F:\Thesis\DNS data\K75_k1';

load 'T.mat';

%%


T=reshape(T,128,128,256);
OH=reshape(OH,128,128,256);
CO=reshape(CO,128,128,256);
CO2=reshape(CO2,128,128,256);
O2=reshape(O2,128,128,256);
H2=reshape(H2,128,128,256);
CH4=reshape(CH4,128,128,256);

%%

clc;

for i=1:128
    for j=1:128
         maxInCO2DNS(i,j)=max(CO2(i,j,:));
         minInCO2DNS(i,j)=min(CO2(i,j,:));
    end
    maxInCO2DNStotal(i)=max(maxInCO2DNS(i,:));
    minInCO2DNStotal(i)=min(minInCO2DNS(i,:));

end

maxCO2BinDNS=max(maxInCO2DNStotal(:));
minCO2BinDNS=min(minInCO2DNStotal(:));
if minCO2BinDNS<0
    minCO2BinDNS=0;
end
clear -regxp maxInCO2DNStotal minInCO2DNStotal maxInCO2DNS minInCO2DNS i j;


%%



numbin=100;
binsize=(maxCO2Bin-minCO2Bin)/numbin;
sumOHDNS=zeros(numbin,1);
checkOHDNS=zeros(numbin,1);
binOHDNS=zeros(numbin,1);
sumCODNS=zeros(numbin,1);
checkCODNS=zeros(numbin,1);
binCODNS=zeros(numbin,1);
sumCO2DNS=zeros(numbin,1);
checkCO2DNS=zeros(numbin,1);
binCO2DNS=zeros(numbin,1);
sumO2DNS=zeros(numbin,1);
checkO2DNS=zeros(numbin,1);
binO2DNS=zeros(numbin,1);
sumH2DNS=zeros(numbin,1);
checkH2DNS=zeros(numbin,1);
binCH4DNS=zeros(numbin,1);
sumCH4DNS=zeros(numbin,1);
checkCH4DNS=zeros(numbin,1);
binH2DNS=zeros(numbin,1);
sumTempDNS=zeros(numbin,1);
checkTempDNS=zeros(numbin,1);
binTempDNS=zeros(numbin,1);
xseiesdataDNS=zeros(numbin,1);

for counter=1:numbin
    F1=(counter-1)*binsize;
    F2=(counter)*binsize;
    
    for i=1:128
        for j=1:128
            for k=1:256
                
                if CO2(i,j,k)>=F1 && CO2(i,j,k)<=F2
                sumOHDNS(counter,1)=sumOHDNS(counter,1)+OH(i,j,k);
                sumCODNS(counter)=sumCODNS(counter,1)+CO(i,j,k);
                sumCO2DNS(counter,1)=sumCO2DNS(counter,1)+CO2(i,j,k);
                sumO2DNS(counter,1)=sumO2DNS(counter,1)+O2(i,j,k);
                sumH2DNS(counter,1)=sumH2DNS(counter,1)+H2(i,j,k);
                sumCH4DNS(counter,1)=sumCH4DNS(counter,1)+CH4(i,j,k);
                sumTempDNS(counter,1)=sumTempDNS(counter,1)+T(i,j,k);

%                 xseiesdataDNS(counter)=((F1+F2)/2)/maxCO2BinDNS;
                
                checkOHDNS(counter,1)=checkOHDNS(counter,1)+1.0;
                checkCODNS(counter,1)=checkCODNS(counter,1)+1.0;
                checkCO2DNS(counter,1)=checkCO2DNS(counter,1)+1.0;
                checkO2DNS(counter,1)=checkO2DNS(counter,1)+1.0;
                checkH2DNS(counter,1)=checkH2DNS(counter,1)+1.0;
                checkCH4DNS(counter,1)=checkCH4DNS(counter,1)+1.0;
                checkTempDNS(counter,1)=checkTempDNS(counter,1)+1.0;
                    
                end
                
            end
        end
    end
    binOHDNS(counter,1)=sumOHDNS(counter,1)/checkOHDNS(counter,1);
    binCODNS(counter,1)=sumCODNS(counter,1)/checkCODNS(counter,1);
    binCO2DNS(counter,1)=sumCO2DNS(counter,1)/checkCO2DNS(counter,1);
    binO2DNS(counter,1)=sumO2DNS(counter,1)/checkO2DNS(counter,1);
    binH2DNS(counter,1)=sumH2DNS(counter,1)/checkH2DNS(counter,1);
    binCH4DNS(counter,1)=sumCH4DNS(counter,1)/checkCH4DNS(counter,1);
    binTempDNS(counter,1)=sumTempDNS(counter,1)/checkTempDNS(counter,1);

end


%%

clear -regxp binOHDNS sumOHDNS checkOHDNS binCODNS sumCODNS checkCODNS; 
clear -regxp binCO2DNS sumCO2DNS checkCO2DNS;
clear -regxp binO2DNS sumO2DNS checkO2DNS binH2DNS sumH2DNS checkH2DNS; 
clear -regxp binCH4DNS sumCH4DNS checkCH4DNS;
clear -regxp binTempDNS sumTempDNS checkTempDNS binsize; 
clear -regxp numbin counter F1 F2 i j;







