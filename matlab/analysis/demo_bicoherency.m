% =============================================================
%            Followed Method Description
% -------------------------------------------------------------


%demo_bicoherence :
%
%Finding the average bicoherence using partitioning in time space
%for each epoch and taking the average  along the ensemble epochs of a
%patient for a specific sleep stage in time space too. Then, we calculate the bicoherence 
%on these averaged values. Finishing by taking the average bicoherence
%for the entire set of patients.
%

% =============================================================
%             Reset your workspace variables
% -------------------------------------------------------------
clear all;

clc;


% =============================================================
%                      Script parameters
% -------------------------------------------------------------



%frequency sampling
fs=256;

%segment duration for partitioning (dt < 30)
dt = 30 ;

%samples per segment (epoch partitioning into K segments)
M = floor(dt * fs) ;

%time vector
t = (0:1:M-1) * (1/fs) ;

%the selected channel 
channel = 'EEG_C4_M1' ;

%length of fft samples 
nfft = 1024;



%the valid sleep stages for patients
sleepstages = {"Sleep stage W";
"Sleep stage N1";
"Sleep stage N2";
"Sleep stage N3";
"Sleep stage R"};


%first patient
p1 = 5 ;

%number of patients
a=10;




% =============================================================
%                     Bicoherence Calculation 
% -------------------------------------------------------------



for j = 1:5 
    
    bicoherency = zeros(nfft,nfft) ; 
    
    
    
    for i=p1:1:p1+a
        
        
        
         S=read_data_of_patient(i);

         
         
        [bic , waxis] = avg_bic(S, channel , sleepstages{j} ,dt , nfft);
        
        
        bicoherency = bicoherency + bic ;
        
        
        
    end
    
    
    bicoherency = bicoherency/(a+1) ;
    
    
    
   
    
   
%    
%     figure(j)
%     clf;
%     contour(waxis,waxis,bicoherency,4), grid on
%     title('Bicoherence estimated via the direct (FFT) method')
%     xlabel('f1'), ylabel('f2')
%   



    %Plot the results
    
    figure(2*(j-1)+1)
    clf
    imagesc(waxis,waxis,bicoherency)
    title(sprintf('Bicoherence estimated via the direct (FFT) method of %s',sleepstages{j}) )
    xlabel('f1'), ylabel('f2')
    colorbar

    
    figure(2*(j-1)+2)
    clf;
    hold on;
    xlabel( 'f1 ');
    ylabel('f2');
    title(sprintf('Bicoherence of  %s', sleepstages{j}));
    s = mesh( -nfft/2:nfft/2-1, -nfft/2:nfft/2-1 , abs(bicoherency(:,:)));
    s.FaceColor = 'interp' ;
    colorbar;
    hold off;

    
    
end