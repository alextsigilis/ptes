%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Function description:
%                 
%--------------------------------------------------------------------------
%This function computes real cepstrum, complex
%cepstrum , complex cepstrum using fft(phase unwrapping  method) , minimum phase 
%impulse response and mixed-phase impulse response (from FFT approach). 
%For these calculations implies signal vector partitioning according to the 
%time length dt. Then, for a specific channel and a selected sleep stage, 
%it partitions epochs into K segments. After, it computes the products for each segment 
%separately  and derives the average products for each epoch record 
%from these computations.In the end , it computes the average products throughout 
%the complete collection of epochs' calculated products.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%
%Function Arguments:
%
%==========================================================================
%Z : the data table of a patient which contains the polysomnographic
%recordings
%channel: the selected channel(EEG_F4_M1 , EEG_C4_M1 , EEG_O2_M1 , EEG_C3_M2 , EMG_chin , EOG_E1_M2 , EOG_E2_M2 , ECG)
%sleep stage: the selected sleep stage()
%zlength : the number of rows of Z
%dt: dt (in seconds) is an integer < 30 seconds  
%nfft: The samples for FFT computation
%--------------------------------------------------------------------------
%
%Return Products:
%
%==========================================================================
%A1:the average computed real cepstrum
%A2:the average computed complex cepstrum
%A3:the average computed complex cepstrum via FFT method (using bicepstral analysis)
%hest_coef:the average estimated mixed-phase impulse response
%hmest_coef:the average estimated minimum-phase impulse response


function [A1, A2, A3, hest_coef ,   hmest_coef] = stageceps(Z, channel, sleep_stage, zlength ,dt,nfft)


%vectors of the function results
rceps_coef = [] ;
cepst_coef1 = [] ;
cepst_coef2 = [] ;
hest_coef = [] ;
hmest_coef = [] ;




% b=hanning(M);


% number of segments per epoch (epoch partitioning)
K = floor(30/dt) ; 





%frequency of sampling
fs=256;

%number of samples for cepstrum estimation using bicepstrum (FFT method)
%algorithm
%nfft = 256;



%samples per segment (epoch partitioning into K segments)
M = floor(dt * fs) ;


c = zeros(M,1) ; 
c1 = zeros(M,1);
c2 = zeros(nfft,1);
h = zeros(nfft,1);

%time vector corresponding to the samples from each segment
t = (0:1:M-1) * (1/fs) ;

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
        
        
       
        
        
        %Applying hanning window on the signal with the aim of smoothing
%         w = hanning(L , 'periodic');
%         
%         A = A .* w ;
        
        
        c = zeros(M,1) ; 
        c1 = zeros(M,1);
        hm = zeros(M,1) ;
        c2 = zeros(nfft,1);
        h = zeros(nfft,1);
        %Apply partitioning for each epoch into K segments and estimate the
        %real and complex cepstrum for each Of them. Then we take the
        %average values of these estimations.
        
            for k = 1:1:K
            
            
            
                 if k*M <= L
                
                

                
                   A11 = A((k-1)*M+1 : k*M);
                   
                   
                   
                   A11 = nonzeros(A11);
                   
                   L1 = length(A11);
                
                   if L1 < 7680 
                      
                       printf("the length of A11 is %d" , L1)
                       
                   end
                   
                   %Estimation of real cepstrum coefficients
                   [cc , hm1] = rceps1(A11);
        
                   c = c + cc ;
        
                   %Estimation of complex cepstrum coefficients using cceps() function of matlab
                   %library
                   c11 = cceps(A11) ;
        
                   c1 = c1 + c11 ;
                   
                   hm = hm + hm1 ;
               
               
                   %Estimation of complex cepstrum coefficients using
                   %bicepstrum (FFT method)
                   [h1,c22,~,~]= bicepsf(A11,20,L1, 0,'unbiased', nfft, 0);
        
                   c2 = c2 + c22 ;
                 
                   h = h + h1 ;
                
                  end
        
%         figure(2)
%         clf;
%         plot(L,A1)
%         
%         
%         A11 = A1-mean(A1);
%         
%         
%         
%        
%         Y=fft(A11);
%         P2 = abs(Y/length(A1))';
%         P1 = P2(1:round(length(A1)/2+1));
%         P1(2:end-1) = 2*P1(2:end-1);
%         
%         
%         
%         figure(3)
%         clf;
%         plot(L,A11)
%         
%         
%         figure(4)
%         clf;
%         hold on;
%         plot(1:length(P1),P1) 
%         title("Single-Sided Amplitude Spectrum of A11(t)")
%         xlabel("SAMPLES")
%         ylabel("|P1(sample)|")
%         hold off;
%         
%         a=max(A11);
%         
%         [p,peaks]=findpeaks(A11,'MinPeakProminence',max(A11)/4);
%         m=mean(diff(peaks));
%         
%         
%         
       
        
           end
        
        
        c  = c./K ;
        c1 = c1./K ;
        c2 = c2./K ;
        h = h./K ;
        hm = hm./K ;
        
        %Store the above estimations in every iteration
        ceps1 = num2cell(c1,1);
        
        ceps2 = num2cell(c2,1);
        
        hest = num2cell(h,1);
        
        rc = num2cell(c,1) ;
        
        hmest = num2cell(hm,1) ;
        
        rceps_coef = [rceps_coef ; rc] ;
        
        cepst_coef1 = [cepst_coef1 ; ceps1];
        
        cepst_coef2 = [cepst_coef2 ; ceps2];
        
        hest_coef = [hest_coef ; hest]  ;
        
        hmest_coef = [hmest_coef ; hmest] ;
        %plot the cepstral coefficients
%         figure(1)
%         clf;
%         subplot(3,1,1)
%         plot(t(1:round(L/2)) , c1(1:round(L/2)) , 'r'); title(sprintf(' Complex Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
% 
% 
%         subplot(3,1,2)
%         plot(t(1:round(L/2)) , c(1:round(L/2)) , 'r'); title(sprintf(' Real Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
%         
% 
%         
%         subplot(3,1,3)
%         plot(-nfft/2:nfft/2-1, c2),title(sprintf(' Complex Cebtrum using bicepsf of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Cepstrum'); grid on
   
           
        
    end
    
end


%Compute the average cepstrum coefficients for every case (real cepstrum , complex cepstrum , complex cepstrum existing from bicepstrum )

A1 = find_avg(rceps_coef);

A2 = find_avg(cepst_coef1);

A3 = find_avg(cepst_coef2);

hest_coef = find_avg(hest_coef);

hmest_coef = find_avg(hmest_coef);

%Plot the average cepstrum coefficients for every case (real cepstrum , complex cepstrum , complex cepstrum existing from bicepstrum )
L = length(A1) ; 

%  figure(2)
%  clf;
%  subplot(4,1,1)
%  plot(t(1:round(L/2)) , A1(1:round(L/2)) , 'r'); title(sprintf('Average Real Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
% 
% 
%  subplot(4,1,2)
%  plot(t(1:round(L/2)) , A2(1:round(L/2)) , 'r'); title(sprintf('Average Complex Cepstrum of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('Time(s)'); ylabel('Cepstrum');
%         
%  subplot(4,1,3)
%  plot(-nfft/2:nfft/2-1, A3),title(sprintf(' Complex Cepstrum using bicepstrum (FFT method) of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Cepstrum'); grid on
% 
%  subplot(4,1,4)
%  plot(-nfft/2:nfft/2-1, hest_coef),title(sprintf(' Impulse response of channel "%s" for the sleep stage %s ', channel , sleep_stage)); xlabel('samples'); ylabel('Impulse response'); grid on
 
