

% Function Description:
% This function can be used to read polysomnographic recordings and
% sleep scorings from EDF files and convert them into
% MATLAB timetables.

% Function Arguments:
% patient: the number of the patient whose recordings
% we want to read. There are 154 patients in the dataset.
% Every patient is given an index ranging from 1 to 154.

% Return Variables:
% S: A timetable which contains the recordings
% and the sleep stage labels for the selected patient.
% This is what the timetable looks like:
%
%  | EEGF4_M1 | EEGC4_M1 | EEGO2_M1 | EEGC3_M2 | EMGChin  | EOGE1_M2 | EOGE2_M2 |   ECG    | Annotations |
% -------------------------------------------------------------------------------------------------------------
%  | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | "Sleep stage W"
%  | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | "Sleep stage W"
%  | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | {7680x1} | "Sleep stage N1"
%
% 
% The following eight columns contain the EEG/EMG/EOG/ECG
% recordings. Every element of this column is a single cell
% containing 30sec x 256Hz = 7680 samples of recorded data.
% The "Annotations" column contains the sleep stage labels
% for every segment. The "Annotations" column contains the sleep stage labels
% for every segment.


function [S] = read_data_of_patient(patient)


if   patient>0 & patient<10

patient_data = sprintf("SN00%d.edf",patient);
patient_annot = sprintf("SN00%d_sleepscoring.edf",patient);

elseif patient<100 
    
    
patient_data = sprintf("SN0%d.edf",patient);
patient_annot = sprintf("SN0%d_sleepscoring.edf",patient);

    
else
    
patient_data = sprintf("SN%d.edf",patient);
patient_annot = sprintf("SN%d_sleepscoring.edf",patient);

    
end



[X,k] = ReadEDF(patient_data);
[r,y] = ReadEDF(patient_annot);



A= string(y.annotation.event)';


B= y.annotation.starttime;

% sampling rate : samples per second
fs = 256;


channels = k.channels;

% segments: the number of records of 1 second duration
segments = k.records;   


channel_names = k.labels ;





valid_stages = cellstr([ ...
    "Sleep stage W",     ...
    "Sleep stage N1",    ...
    "Sleep stage N2",    ...
    "Sleep stage N3",    ...
    "Sleep stage R"]);


rows = ~ismember(                ...
    categorical(A),  ...
    categorical(valid_stages)    ...
    );

% duration in seconds of every segment
w = 30;

% samples per segment
samp = w * fs ; 


%remove rows from A,B that don't belong to valid stages of sleep
A(rows) = [];
B(rows) = [];





S = cell(length(A),channels);



for i=1:1:8
    
  for j=1:1:length(B)
      
      
      x1 = X{1,i}((samp*(j-1)+1):samp*j,1);
      x1 = num2cell(x1,1);
      S{j,i} = x1{1,1} ;
      
  end
  
end

    
    


S = cell2table(S);
S.Properties.VariableNames{1,1} = 'EEG_F4_M1';
S.Properties.VariableNames{1,2} = 'EEG_C4_M1';
S.Properties.VariableNames{1,3} = 'EEG_O2_M1';
S.Properties.VariableNames{1,4} = 'EEG_C3_M2';
S.Properties.VariableNames{1,5} = 'EMG_chin';
S.Properties.VariableNames{1,6} = 'EOG_E1_M2' ; 
S.Properties.VariableNames{1,7} = 'EOG_E2_M2' ;
S.Properties.VariableNames{1,8} = 'ECG' ;

S = addvars(S,A,'NewVariableNames','Annotations');
S = addvars(S,B,'NewVariableNames','Start_time');
S = addvars(S,B+30,'NewVariableNames','End_time');
