
% =============================================================
%            Followed Method Description
% -------------------------------------------------------------

%
%
%demo_cepstrum2 :
%
%Finding the average cepstral coefficients using partitioning in time space
%for each epoch and  calculating the cepstral coefficients for these
%computed values. Then, we take the average of calculated cepstral
%coefficients for the ensemble epochs. Finishing by taking the average cepstral
%coefficients for the entire set of patients.
%




% =============================================================
%             Reset your workspace variables
% -------------------------------------------------------------
clear all;

clc;

% =============================================================
%                      Save the results 
% -------------------------------------------------------------

%define the name of user  
USER = 'giorg';

% save directory
dir = sprintf("C:\\Users\\%s\\Desktop",USER);

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\Cepstrum_statistics1" ,dir               ...
        ));
    
mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\Cepstrum_graphs1" ,dir               ...
        ));




% =============================================================
%                      Script parameters
% -------------------------------------------------------------


%frequency sampling
fs=256;

%segment duration for partitioning (dt < 30)
dt = 30 ;

%samples per segment (epoch partitioning into K segments , K=30/dt)
M = floor(dt * fs) ;


%time vector
t = (0:1:M-1) * (1/fs) ;

%selected channel from the list: 'EEG_F4_M1' , 'EEG_C4_M1' , 'EEG_O2_M1' ,
%'EEG_C3_M2' , 'EMG_chin' ,'EOG_E1_M2' , 'EOG_E2_M2' , 'ECG'
channel = {'EEG_F4_M1' , 'EEG_C4_M1' , 'EEG_O2_M1' ,'EEG_C3_M2' , 'EMG_chin' ,'EOG_E1_M2' , 'EOG_E2_M2' , 'ECG'} ;

%the fft samples
nfft = 256;

%***!Select the p1 (patient 1) and a from the valid patient's data:
%[1:13,15:17,19:31,37:63,65:76,78:97,99:111,113:119]!*****

%number of patients
a=10;                      %a=0 for using the whole set of patient's data 

%first patient
p1 = 1 ;

%table for storing cepstral coefficients and impulse response estimations
names = ["real_cepstrum_coef" "complex_ceps_coeff" "complex_ceps_coeff_fft_method" "minimum_phase_impulse_response" "mixed_phase_impulse_response_fft_method" "Patient" "channel" "Annotations"];
types = ["cell" "cell" "cell" "cell" "cell" "double" "string" "string"];
sz = [5*(a+1) numel(types)];
cepstrum_coeff = table('Size',sz,'VariableTypes',types,'VariableNames',names);

%table for storing statistics
names = ["var" "skw" "krt" "entropia" "Patient" "channel" "Annotations"];
types = [ "double" "double" "double" "double" "double" "string" "string"];
sz = [5*(a+1) numel(types)];
cepstrum_statistics = table('Size',sz,'VariableTypes',types,'VariableNames',names);

%the valid sleep stages for patients
sleepstages = {"Sleep stage W";
"Sleep stage N1";
"Sleep stage N2";
"Sleep stage N3";
"Sleep stage R"};





% =============================================================
%                      Cepstrum Calculations
% -------------------------------------------------------------

for k=1:numel(channel)
    
    
    
    for j = 1:5 
    
        %real cepstrum coefficients
        real_cepstrum = zeros(M,1) ; 
    
        %complex cepstrum coefficients
        complex_cepstrum1 = zeros(M,1) ;
    
        %estimated minimum phase impulse response
        hm = zeros(M,1) ; 
    
        %complex cepstrum coefficients using fft method
        complex_cepstrum2 = zeros(1,nfft) ;
    
        %estimated mixed-phase impulse response using fft method
        h = zeros(1,nfft) ;
    
    
        %a=0;
        
        %Find the average of the cepstral calculated coefficients for the entire set of patients
    
        for i=p1:p1+a-1
        
        
           %data of patient
           S=read_data_of_patient(i);

         
           %cepstral and impulse response estimated coefficients of patient 
           [RC,CC1,CC2 , H , Hm] = stageceps(S, channel{k} , sleepstages{j} , height(S),dt,nfft);
        
           %storing cepstral coefficients and impulse response estimations
           cepstrum_coeff(5*(118+1)*(k-1)+(i-p1)*5+j,:) = {RC CC1 CC2 H Hm i channel{k} sleepstages{j}};
        
        
           %computing cepstrum statistics
           cepstrum_statistics(5*(118+1)*(k-1)+(i-p1)*5+j,:) = {std(CC1) skewness(CC1) kurtosis(CC1)  abs(mean(CC1 .* log2(CC1),'omitnan')) i channel{k} sleepstages{j}};
        
        
           real_cepstrum = real_cepstrum + RC ;
        
           complex_cepstrum1 = complex_cepstrum1 + CC1 ;
        
           complex_cepstrum2 = complex_cepstrum2 + CC2 ;
        
           h = h + H ;
        
           hm = hm + Hm ;
           
           
           %Because the patient's data are not accesible in a continuous
           %manner, a change is made for every iteration of counting the total
           %number of patients
           %a=a+1;
        
        
        end
    
       %find the average coefficients for the set of patients
    
       real_cepstrum = real_cepstrum/(a) ;
    
       complex_cepstrum1 = complex_cepstrum1/(a) ;
    
       complex_cepstrum2 = complex_cepstrum2/(a)  ;
    
       h = h/(a)  ;
    
       hm = hm/(a) ;
    
   
    
       real_cepstrum = fftshift(real_cepstrum);
    
       complex_cepstrum1 = fftshift(complex_cepstrum1);
    
   
   
       %Plot the results
    
       figure(5*(k-1)+j)
       clf;
       subplot(5,1,1)
       plot((-M/2:(M/2-1))*(1/fs) , real_cepstrum(1:round(M)) , 'r'); title(sprintf('Average Real Cepstrum of channel "%s" for the %s ', channel{k} , sleepstages{j})); xlabel('Time(s)'); ylabel('Real Cepstrum');


       subplot(5,1,2)
       plot((-M/2:(M/2-1))*(1/fs) , complex_cepstrum1(1:round(M)) , 'r'); title(sprintf('Average Complex Cepstrum of channel "%s" for the %s ', channel{k} , sleepstages{j})); xlabel('Time(s)'); ylabel('Complex Cepstrum');
     
       subplot(5,1,3)
       plot(t(1:round(M/2)) , hm(1:round(M/2)) , 'r'); title(sprintf('Minimum phase impulse response of channel "%s" for the %s ', channel{k} , sleepstages{j})); xlabel('Time(s)'); ylabel('Minimum phase impulse response');
    
    
       subplot(5,1,4)
       plot(-nfft/2:nfft/2-1, complex_cepstrum2),title(sprintf(' Complex Cepstrum (FFT method) of channel "%s" for the %s ', channel{k} , sleepstages{j})); xlabel('samples'); ylabel('Cepstrum'); grid on

       subplot(5,1,5)
       plot(-nfft/2:nfft/2-1, h),title(sprintf(' Impulse response of channel "%s" for the %s ', channel{k} , sleepstages{j})); xlabel('samples'); ylabel('Impulse response'); grid on
    
       
       filename = sprintf("%s\\Cepstrum_graphs1\\%d.fig",dir, 5*(k-1)+2*j);
       saveas(gcf, filename);
    
    end
    
end


filename = sprintf("%s\\Cepstrum_statistics1\\cepstrum_statistics_of_patients1.mat", dir);
save(filename,"cepstrum_statistics")


filename = sprintf("%s\\Cepstrum_statistics1\\cepstrum_coefficients_of_patients1.mat", dir);
save(filename,"cepstrum_coeff")


% =============================================================
%                      Plot Statistics
% -------------------------------------------------------------

%Plotting a histogram of the variance of cepstral coefficients for each slape stage of
%the selected channel

for kk=1:8
    
figure(5*(k-1)+j+kk)
clf

    for jj=1:5
    
        a1=[];
    
        for i=[1:13,15:17,19:35,37:63,65:76,78:97,99:119]
            
            st = double(table2array(cepstrum_statistics(5*(a)*(kk-1)+(i-p1)*5+jj,1)));
            
            if st ~= 0
    
            a1=[a1; st];
            
            end
        
        end
    
        edges = linspace(0, 15 ,round(sqrt(numel(a1))));
    
        hold on;
        histogram(a1,edges);

    end


legend(["W", "N1", "N2", "N3", "R"]);
xlabel("Variance of Cepstrum ")
ylabel("Counts")
title(sprintf("Variance of Cepstrum for channel %s",channel{kk}))
hold off;


filename = sprintf("%s\\Cepstrum_statistics1\\Variance_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);
end 



%Plotting a histogram of the skewness of cepstral coefficients for each slape stage of
%the selected channel

for kk=1:8
    
figure(5*(k-1)+j+kk+8)
clf

    for jj=1:5
    
        a1=[];
    
        for i=[1:13,15:17,19:35,37:63,65:76,78:97,99:119]
            
            st = double(table2array(cepstrum_statistics(5*(a)*(kk-1)+(i-p1)*5+jj,2)));
            
            if st ~= 0
    
            a1=[a1; st];
            
            end
        
        end
    
        edges = linspace(0, 5 ,round(sqrt(numel(a1))));
    
        hold on;
        histogram(a1,edges);

    end


legend(["W", "N1", "N2", "N3", "R"]);
xlabel("Skewness of Cepstrum")
ylabel("Counts")
title(sprintf("Skewness of Cepstrum for channel %s ",channel{kk}))
hold off;

filename = sprintf("%s\\Cepstrum_statistics1\\Skewness_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);

end 


%Plotting a histogram of the kyrtosis of cepstral coefficients for each slape stage of
%the selected channel



for kk=1:8
    
figure(5*(k-1)+j+kk+16)
clf

    for jj=1:5
    
        a1=[];
    
        for i=[1:13,15:17,19:35,37:63,65:76,78:97,99:119]
            
            st = double(table2array(cepstrum_statistics(5*(a)*(kk-1)+(i-p1)*5+jj,3)));
            
            if st ~= 0
    
            a1=[a1; st];
            
            end
        
        end
    
        edges = linspace(0, 4000 ,round(sqrt(numel(a1))));
    
        hold on;
        histogram(a1,edges);

    end


legend(["W", "N1", "N2", "N3", "R"]);
xlabel("Kyrtosis of Cepstrum")
ylabel("Counts")
title(sprintf(" Kyrtosis of Cepstrum for channel %s",channel{kk}))
hold off;

filename = sprintf("%s\\Cepstrum_statistics1\\Kyrtosis_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);
end 



%Plotting a histogram of the entropia of cepstral coefficients for each slape stage of
%the selected channel


for kk=1:8
    
figure(5*(k-1)+j+kk+24)
clf

    for jj=1:5
    
        a1=[];
    
        for i=[1:13,15:17,19:35,37:63,65:76,78:97,99:119]
            
            st = double(table2array(cepstrum_statistics(5*(a)*(kk-1)+(i-p1)*5+jj,4)));
            
            if  st ~= 0
    
            a1=[a1; st];
            
            end
        
        end
    
        edges = linspace(0, 1 ,round(sqrt(numel(a1))));
    
        hold on;
        histogram(a1,edges);

    end


legend(["W", "N1", "N2", "N3", "R"]);
xlabel("Entropia of Cepstrum")
ylabel("Counts")
title(sprintf(" Entropia of Cepstrum for channel %s",channel{kk}))
hold off;

filename = sprintf("%s\\Cepstrum_statistics1\\Entropia_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);

end 