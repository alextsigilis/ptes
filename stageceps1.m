function [A1, A2, A3, hest_coef] = stageceps1(Z, channel, sleep_stage, zlength )


%vectors of the function results
rceps_coef = [] ;
cepst_coef1 = [] ;
cepst_coef2 = [] ;
cepst_coef3 = [] ;
hest_coef = [] ;
hest_coef1 = [] ;


%frequency of sampling
fs=256;

%number of samples for cepstrum estimation using bicepstrum (FFT method)
%algorithm
nfft = 256;


% channel matching with the number of the column in Z matrix
if channel == "EEG_F4_M1"
    
    j = 1;
    
elseif   channel == "EEG_C4_M1"
    
     j = 2;
    
elseif   channel == "EEG_O2_M1"
    
     j = 3;
    
elseif   channel == "EEG_C3_M2"
    
     j = 4;
    
elseif   channel == "EMG_chin"
    
     j = 5;
    
elseif   channel == "EOG_E1_M2"
    
     j = 6;
    
elseif   channel == "EOG_E2_M2"
    
     j = 7;
    
else
    
     j = 8;
    
end


% Find real and complex cepstrum coefficients for a specific sleep stage
% from the channel "j" selected above
% and store them for future processing
for i = 1:zlength
    
   if string(Z.Annotations{i}) == sleep_stage 
        
        
        %A is the epoch segment for the specific sleep stage and channel j
        A = Z{i,j}{1,1};
        
        L=length(A);
        
        
        %time vector corresponding to the samples from each segment
        t = (0:1:L-1) * (1/fs) ;
        
        
        %Applying hanning window on the signal with the aim of smoothing
        w = hanning(L , 'periodic');
        
        A = A .* w ;
        
        %Estimation of real cepstrum coefficients
        [c , hm] = rceps(A);
        
                
        
       %Estimation of complex cepstrum coefficients using cceps() function of matlab
       %library
       c1 = cceps(A) ;
        
                
               
               
       %Estimation of complex cepstrum coefficients using
       %bicepstrum (FFT method)
       [h,c2,~,~]= bicepsf(A,20,L, 0,'unbiased', nfft, 0);
       
       
       [h1,c3,~,~,~,~] = biceps(A,10,10,L,0,'unbiased',2*(10+10));
       
        %Store the above estimations in every iteration
        ceps1 = num2cell(c1,1);
        
        ceps2 = num2cell(c2,1);
        
        ceps3 = num2cell(c3,1);
        
        hest = num2cell(h,1);
        
        hest1 = num2cell(h1,1);
        
        rc = num2cell(c,1) ;
        
        rceps_coef = [rceps_coef ; rc] ;
        
        cepst_coef1 = [cepst_coef1 ; ceps1];
        
        cepst_coef2 = [cepst_coef2 ; ceps2];
        
        cepst_coef3 = [cepst_coef3 ; ceps3];
        
        hest_coef = [hest_coef ; hest]  ;
        
        hest_coef1 = [hest_coef1 ; hest1]  ;
    end
end

%Compute the average cepstrum coefficients for every case (real cepstrum , complex cepstrum , complex cepstrum existing from bicepstrum )

A1 = find_avg(rceps_coef);

A2 = find_avg(cepst_coef1);

A3 = find_avg(cepst_coef2);

A4 = find_avg(cepst_coef3);

hest_coef = find_avg(hest_coef);

hest_coef1 = find_avg(hest_coef1);


%Plot the average cepstrum coefficients for every case (real cepstrum , complex cepstrum , complex cepstrum existing from bicepstrum )
L = length(A1) ; 

 figure(2)
 clf;
 subplot(4,1,1)
 plot(t(1:round(L/2)) , A1(1:round(L/2)) , 'r'); title(sprintf('Average Complex Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');


 subplot(4,1,2)
 plot(t(1:round(L/2)) , A2(1:round(L/2)) , 'r'); title(sprintf('Average Real Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
        
 subplot(4,1,3)
 plot(-nfft/2:nfft/2-1, A3),title(sprintf(' Complex Cepstrum using bicepstrum (FFT method) of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Cepstrum'); grid on

 subplot(4,1,4)
 plot(-nfft/2:nfft/2-1, hest_coef),title(sprintf(' Impulse response of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Impulse response'); grid on
 
figure(3)
clf;
subplot(2,1,1)
plot(1:length(A4) , A4 , 'r'); title(sprintf('Average Complex Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
        
       
subplot(2,1,2)       
plot(1:length(hest_coef1), hest_coef1),title(sprintf(' Impulse response of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Impulse response'); grid on 
 
        