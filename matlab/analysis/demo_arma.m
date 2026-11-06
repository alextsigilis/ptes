% =============================================================
%            Followed Method Description
% -------------------------------------------------------------


%demo_arma :
%
%Finding the average estimated arma coefficients using partitioning in time space
%for each epoch and taking the average  along the ensemble epochs of a
%patient for a specific sleep stage in time space too. Then, we estimate the orders p ,q 
%of ar and ma part of arma model respectively on these averaged values.
%Storing the arma model computed coefficients and
%their orders in a table for the selected channel of each patient. 
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

dt = 30 ;

%samples per segment (epoch partitioning into K segments)
M = floor(dt * fs) ;

%time vector
t = (0:1:M-1) * (1/fs) ;

%selected channel from the list: 'EEG_F4_M1' , 'EEG_C4_M1' , 'EEG_O2_M1' ,
%'EEG_C3_M2' , 'EMG_chin' ,'EOG_E1_M2' , 'EOG_E2_M2' , 'ECG' 
channel = 'EOG_E1_M2' ;


% channels = {'EEG_F4_M1';
% 'EEG_C4_M1';
% 'EEG_O2_M1';
% 'EEG_C3_M2';
% 'EMG_chin';
% 'EOG_E1_M2';
% 'EOG_E2_M2';
% 'ECG'};

%the valid sleep stages for patients
sleepstages = {"Sleep stage W";
"Sleep stage N1";
"Sleep stage N2";
"Sleep stage N3";
"Sleep stage R"};

%number of patients
a=10;



 names = ["Sleep_stage_W","Sleep_stage_N1","Sleep_stage_N2","Sleep_stage_N3","Sleep_stage_R","Patient","channel"];
 types = ["cell","cell","cell","cell","cell","double","string"];
 sz = [a+1 numel(types)];
 armamodel = table('Size',sz,'VariableTypes',types,'VariableNames',names);
 

%first patient
p1 = 5 ;



% =============================================================
%                     Arma model Calculation 
% -------------------------------------------------------------


for j = 1:5 
    
    
    
    
    
    for i=p1:1:p1+a
        
        
         %patient data vector
         S=read_data_of_patient(i);

         
        %ar and ma estimated coefficients and theirs orders p , q respectively 
        [ar , ma ,p , q] = arma_estimate(S,channel, sleepstages{j},dt);
        
        %Store the computed values
        C = {{ar , ma , p ,q}} ;
        
        
        
        armamodel{i-p1+1,j} = C ;
        armamodel.Patient(i-p1+1) = i ;
        armamodel.channel(i-p1+1) = channel ;
        
        
    end
    
       

    
end