% =============================================================
%            Followed Method Description
% -------------------------------------------------------------


%demo_bicepstrum :
%
%Finding the average bicepstral coefficients using partitioning in time space
%for each epoch and taking the average  along the ensemble epochs of a
%patient for a specific sleep stage in time space too. Then, we calculate
%the bicepstral coefficients on these averaged values. Finishing by taking
%the average bicepstral coefficients for the entire set of patients.
%

% =============================================================
%             Reset your workspace variables
% -------------------------------------------------------------
clear all;

clc;


% =============================================================
%                      Save Results
% -------------------------------------------------------------

%set the USER 
USER = 'giorg';

% save directory
dir = sprintf("C:\\Users\\%s\\Desktop",USER);

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\Bicepstrum_statistics" ,dir               ...
        ));
    
mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\Bicepstrum_graphs" ,dir               ...
        ));    
    
 
    
% =============================================================
%                      Script parameters
% -------------------------------------------------------------
    
    
%frequency sampling per epoch
fs=256;

%fft samples
nfft = 256;

%frequency vector by taking 
om = [-nfft/2:nfft/2-1] / nfft;

%selected channel from the list: 'EEG_F4_M1' , 'EEG_C4_M1' , 'EEG_O2_M1' ,
%'EEG_C3_M2' , 'EMG_chin' ,'EOG_E1_M2' , 'EOG_E2_M2' , 'ECG'
channel = {'EEG_F4_M1' , 'EEG_C4_M1' , 'EEG_O2_M1' ,'EEG_C3_M2' , 'EMG_chin' ,'EOG_E1_M2' , 'EOG_E2_M2' , 'ECG' }; 



B = cell(5,1);

m=zeros(1,5);

%the valid sleep stages
sleepstages = {"Sleep stage W";
"Sleep stage N1";
"Sleep stage N2";
"Sleep stage N3";
"Sleep stage R"};


%first patient
p1 = 15 ;

%last patient
p2=15;

%number of patients : a+1
a=0;

%table for storing statistics
names = ["Bicepstrum" "var" "skw" "krt" "entropia" "squared_entropy" "moment2nd" "Patient"  "channel"  "Annotations"];
types = [ "cell" "cell" "cell" "cell" "cell" "cell" "cell" "double" "string" "string"];
sz = [5*(p2-p1+1) numel(types)];
Bicepstrum_statistics = table('Size',sz,'VariableTypes',types,'VariableNames',names);








% =============================================================
%                      Bicepstrum Calculation
% -------------------------------------------------------------



for k=1:8 
    
    
    
  for j = 1:5 
    
      Bcepsmat = zeros(nfft) ;
      
      a=0;
      

         for i=p1:p2 %i=[1:13,15:17,19:31,37:63,65:76,78:97,99:111,113:119] 





             S=read_data_of_patient(i);

             A = S{3,1}{1,1} ;



            %[B] = Bicepstrum(A, 20 , length(A) , 41) ;


            %calculating bicepstral coefficients
            [Bceps,m1]=avg_bicepstrum(S,channel{k}, sleepstages{j},nfft,30);
        
            %computing and storing Bicepstrum statistics
            Bicepstrum_statistics(5*(p2-p1+1)*(k-1)+(i-p1)*5+j,:) = {Bceps std(Bceps)' skewness(Bceps)' kurtosis(Bceps)'  abs(mean(Bceps .* log2(Bceps),'omitnan'))' abs(mean(sqrt(Bceps) .* log2(sqrt(Bceps)),'omitnan'))' m1 i channel{k} sleepstages{j}};
 
            Bcepsmat = Bcepsmat + Bceps ;
            
            a=a+1;
        
        

         end

          %computing the average bicepstral coefficients
          Bcepsmat = Bcepsmat/(a) ;
          %B(j) = num2cell(Bcepsmat ,[1 2]);

          figure(10*(k-1)+2*(j-1)+1)
          clf;
          hold on;
          xlabel( 'Lag m ');
          ylabel('Lag n');
          title(sprintf('Bicepstrum of  %s for channel %s', sleepstages{j} , channel{k}));
          s = mesh( -nfft/2:nfft/2-1, -nfft/2:nfft/2-1 , abs(Bcepsmat(:,:)));
          s.FaceColor = 'interp' ;
          colorbar;
          hold off;
          
          filename = sprintf("%s\\Bicepstrum_graphs\\%d.png",dir, 10*(k-1)+2*(j-1)+1);
          saveas(gcf, filename);

          figure(10*(k-1)+2*(j-1)+2)
          clf;
          clf, contour(om,om,abs(Bcepsmat),8), grid,
          title(sprintf('Bicepstrum of  %s for channel %s', sleepstages{j} , channel{k}));
          
          filename = sprintf("%s\\Bicepstrum_graphs\\%d.png",dir, 10*(k-1)+2*(j-1)+2);
          saveas(gcf, filename);

 


%m(j) = mom2(B{j},om);


   end

end

filename = sprintf("%s\\Bicepstrum_statistics\\bicepstrum_stats.mat", dir);
save(filename,"Bicepstrum_statistics")


% =============================================================
%                      Plot Statistics
% -------------------------------------------------------------



%Make histogram of Bicepstrum variance for sleep stages W , R



for kk=1:8
    
figure(10*(k-1)+2*(5-1)+2+kk)
clf

      
    
          
          a1 = cell2mat(Bicepstrum_statistics.var(Bicepstrum_statistics.Annotations == 'Sleep stage W' & Bicepstrum_statistics.channel == channel{kk}));
          a2 = cell2mat(Bicepstrum_statistics.var(Bicepstrum_statistics.Annotations == 'Sleep stage R' & Bicepstrum_statistics.channel == channel{kk}));
          a3 = cell2mat(Bicepstrum_statistics.var(Bicepstrum_statistics.Annotations == 'Sleep stage N1' & Bicepstrum_statistics.channel == channel{kk}));
          a4 = cell2mat(Bicepstrum_statistics.var(Bicepstrum_statistics.Annotations == 'Sleep stage N2' & Bicepstrum_statistics.channel == channel{kk}));
          a5 = cell2mat(Bicepstrum_statistics.var(Bicepstrum_statistics.Annotations == 'Sleep stage N3' & Bicepstrum_statistics.channel == channel{kk}));
%           edges1 = linspace(0, max(a1) ,4*max(a1));
%           edges2 = linspace(0,max(a2),4*max(a2));
%           edges3 = linspace(0,max(a3),4*max(a3));
%           edges4 = linspace(0,max(a4),4*max(a4));
%           edges5 = linspace(0,max(a5),4*max(a5));
    
          hold on;
          histogram(a1);
          histogram(a2);
          histogram(a3);
          histogram(a4);
          histogram(a5);

       


%legend(["W", "N1", "N2", "N3", "R"]);
legend(["W" , "R" , "N1" , "N2" , "N3"]);
xlabel("Variance of Bicepstrum")
ylabel("Counts")
title(sprintf("Variance of |B_3^x(f_1, f_2)| of channel %s",channel{kk}))
hold off;
filename = sprintf("%s\\Bicepstrum_statistics\\Variance_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);


end





%Make histogram of Bicepstrum skewness for sleep stages W , R

for kk=1:8
    
figure(10*(k-1)+2*(5-1)+10+kk)
clf

      
    
          
          a1 = cell2mat(Bicepstrum_statistics.skw(Bicepstrum_statistics.Annotations == 'Sleep stage W' & Bicepstrum_statistics.channel == channel{kk}));
          a2 = cell2mat(Bicepstrum_statistics.skw(Bicepstrum_statistics.Annotations == 'Sleep stage R' & Bicepstrum_statistics.channel == channel{kk}));
          a3 = cell2mat(Bicepstrum_statistics.skw(Bicepstrum_statistics.Annotations == 'Sleep stage N1' & Bicepstrum_statistics.channel == channel{kk}));
          a4 = cell2mat(Bicepstrum_statistics.skw(Bicepstrum_statistics.Annotations == 'Sleep stage N2' & Bicepstrum_statistics.channel == channel{kk}));
          a5 = cell2mat(Bicepstrum_statistics.skw(Bicepstrum_statistics.Annotations == 'Sleep stage N3' & Bicepstrum_statistics.channel == channel{kk}));
%           edges1 = linspace(0, max(a1) ,4*max(a1));
%           edges2 = linspace(0,max(a2),4*max(a2));
%           edges3 = linspace(0,max(a3),4*max(a3));
%           edges4 = linspace(0,max(a4),4*max(a4));
%           edges5 = linspace(0,max(a5),4*max(a5));
    
          hold on;
          histogram(a1);
          histogram(a2);
          histogram(a3);
          histogram(a4);
          histogram(a5);

       


%legend(["W", "N1", "N2", "N3", "R"]);
legend(["W" , "R" , "N1" , "N2" , "N3"]);
xlabel("Skewness of Bicepstrum")
ylabel("Counts")
title(sprintf("Skewness of |B_3^x(f_1, f_2)| of channel %s",channel{kk}))
hold off;
filename = sprintf("%s\\Bicepstrum_statistics\\Skewness_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);


end





%Make histogram of Bicepstrum kyrtosis for sleep stages W , R


for kk=1:8
    
figure(10*(k-1)+2*(5-1)+18+kk)
clf

      
    
          
          a1 = cell2mat(Bicepstrum_statistics.krt(Bicepstrum_statistics.Annotations == 'Sleep stage W' & Bicepstrum_statistics.channel == channel{kk}));
          a2 = cell2mat(Bicepstrum_statistics.krt(Bicepstrum_statistics.Annotations == 'Sleep stage R' & Bicepstrum_statistics.channel == channel{kk}));
          a3 = cell2mat(Bicepstrum_statistics.krt(Bicepstrum_statistics.Annotations == 'Sleep stage N1' & Bicepstrum_statistics.channel == channel{kk}));
          a4 = cell2mat(Bicepstrum_statistics.krt(Bicepstrum_statistics.Annotations == 'Sleep stage N2' & Bicepstrum_statistics.channel == channel{kk}));
          a5 = cell2mat(Bicepstrum_statistics.krt(Bicepstrum_statistics.Annotations == 'Sleep stage N3' & Bicepstrum_statistics.channel == channel{kk}));
%           edges1 = linspace(0, max(a1) ,4*max(a1));
%           edges2 = linspace(0,max(a2),4*max(a2));
%           edges3 = linspace(0,max(a3),4*max(a3));
%           edges4 = linspace(0,max(a4),4*max(a4));
%           edges5 = linspace(0,max(a5),4*max(a5));
    
          hold on;
          histogram(a1);
          histogram(a2);
          histogram(a3);
          histogram(a4);
          histogram(a5);

       


%legend(["W", "N1", "N2", "N3", "R"]);
legend(["W" , "R" , "N1" , "N2" , "N3"]);
xlabel("Kyrtosis of Bicepstrum")
ylabel("Counts")
title(sprintf("Kyrtosis of |B_3^x(f_1, f_2)| of channel %s",channel{kk}))
hold off;
filename = sprintf("%s\\Bicepstrum_statistics\\Kyrtosis_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);


end





%Make histogram of Bicepstrum Entropia for sleep stages W , R


for kk=1:8
    
figure(10*(k-1)+2*(5-1)+26+kk)
clf

      
    
          
          a1 = cell2mat(Bicepstrum_statistics.entropia(Bicepstrum_statistics.Annotations == 'Sleep stage W' & Bicepstrum_statistics.channel == channel{kk}));
          a2 = cell2mat(Bicepstrum_statistics.entropia(Bicepstrum_statistics.Annotations == 'Sleep stage R' & Bicepstrum_statistics.channel == channel{kk}));
          a3 = cell2mat(Bicepstrum_statistics.entropia(Bicepstrum_statistics.Annotations == 'Sleep stage N1' & Bicepstrum_statistics.channel == channel{kk}));
          a4 = cell2mat(Bicepstrum_statistics.entropia(Bicepstrum_statistics.Annotations == 'Sleep stage N2' & Bicepstrum_statistics.channel == channel{kk}));
          a5 = cell2mat(Bicepstrum_statistics.entropia(Bicepstrum_statistics.Annotations == 'Sleep stage N3' & Bicepstrum_statistics.channel == channel{kk}));
%           edges1 = linspace(0, max(a1) ,4*max(a1));
%           edges2 = linspace(0,max(a2),4*max(a2));
%           edges3 = linspace(0,max(a3),4*max(a3));
%           edges4 = linspace(0,max(a4),4*max(a4));
%           edges5 = linspace(0,max(a5),4*max(a5));
    
          hold on;
          histogram(a1);
          histogram(a2);
          histogram(a3);
          histogram(a4);
          histogram(a5);

       


%legend(["W", "N1", "N2", "N3", "R"]);
legend(["W" , "R" , "N1" , "N2" , "N3"]);
xlabel("Entropia of Bicepstrum")
ylabel("Counts")
title(sprintf("Entropia of |B_3^x(f_1, f_2)| of channel %s",channel{kk}))
hold off;
filename = sprintf("%s\\Bicepstrum_statistics\\Entropia_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);


end


for kk=1:8
    
figure(10*(k-1)+2*(5-1)+34+kk)
clf

      
    
          
          a1 = cell2mat(Bicepstrum_statistics.squared_entropy(Bicepstrum_statistics.Annotations == 'Sleep stage W' & Bicepstrum_statistics.channel == channel{kk}));
          a2 = cell2mat(Bicepstrum_statistics.squared_entropy(Bicepstrum_statistics.Annotations == 'Sleep stage R' & Bicepstrum_statistics.channel == channel{kk}));
          a3 = cell2mat(Bicepstrum_statistics.squared_entropy(Bicepstrum_statistics.Annotations == 'Sleep stage N1' & Bicepstrum_statistics.channel == channel{kk}));
          a4 = cell2mat(Bicepstrum_statistics.squared_entropy(Bicepstrum_statistics.Annotations == 'Sleep stage N2' & Bicepstrum_statistics.channel == channel{kk}));
          a5 = cell2mat(Bicepstrum_statistics.squared_entropy(Bicepstrum_statistics.Annotations == 'Sleep stage N3' & Bicepstrum_statistics.channel == channel{kk}));
%           edges1 = linspace(0, max(a1) ,4*max(a1));
%           edges2 = linspace(0,max(a2),4*max(a2));
%           edges3 = linspace(0,max(a3),4*max(a3));
%           edges4 = linspace(0,max(a4),4*max(a4));
%           edges5 = linspace(0,max(a5),4*max(a5));
    
          hold on;
          histogram(a1);
          histogram(a2);
          histogram(a3);
          histogram(a4);
          histogram(a5);

       


%legend(["W", "N1", "N2", "N3", "R"]);
legend(["W" , "R" , "N1" , "N2" , "N3"]);
xlabel("Squared Entropy of Bicepstrum")
ylabel("Counts")
title(sprintf("Squared Entropy of |B_3^x(f_1, f_2)| of channel %s",channel{kk}))
hold off;
filename = sprintf("%s\\Bicepstrum_statistics\\Squared_Entropia_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);


end



%Make histogram of Bicepstrum 2nd_moment for sleep stages W , R

for kk=1:8
    
figure(10*(k-1)+2*(5-1)+42+kk)
clf

      
    
          
          a1 = cell2mat(Bicepstrum_statistics.moment2nd(Bicepstrum_statistics.Annotations == 'Sleep stage W' & Bicepstrum_statistics.channel == channel{kk}));
          a2 = cell2mat(Bicepstrum_statistics.moment2nd(Bicepstrum_statistics.Annotations == 'Sleep stage R' & Bicepstrum_statistics.channel == channel{kk}));
          a3 = cell2mat(Bicepstrum_statistics.moment2nd(Bicepstrum_statistics.Annotations == 'Sleep stage N1' & Bicepstrum_statistics.channel == channel{kk}));
          a4 = cell2mat(Bicepstrum_statistics.moment2nd(Bicepstrum_statistics.Annotations == 'Sleep stage N2' & Bicepstrum_statistics.channel == channel{kk}));
          a5 = cell2mat(Bicepstrum_statistics.moment2nd(Bicepstrum_statistics.Annotations == 'Sleep stage N3' & Bicepstrum_statistics.channel == channel{kk}));
%           edges1 = linspace(0, max(a1) ,4*max(a1));
%           edges2 = linspace(0,max(a2),4*max(a2));
%           edges3 = linspace(0,max(a3),4*max(a3));
%           edges4 = linspace(0,max(a4),4*max(a4));
%           edges5 = linspace(0,max(a5),4*max(a5));
    
          hold on;
          histogram(a1);
          histogram(a2);
          histogram(a3);
          histogram(a4);
          histogram(a5);

       


%legend(["W", "N1", "N2", "N3", "R"]);
legend(["W" , "R" , "N1" , "N2" , "N3"]);
xlabel("2nd_moment of Bicepstrum")
ylabel("Counts")
title(sprintf("2nd_moment of |B_3^x(f_1, f_2)| of channel %s",channel{kk}))
hold off;
filename = sprintf("%s\\Bicepstrum_statistics\\2nd_moment_histogram_for_%s.png", dir, channel{kk});
saveas(gcf, filename);


end