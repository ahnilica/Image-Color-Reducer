% ECE471/571 project 1
% How to illustrate data points
%

% clear the figure
clf;

% load the training set
load RESULT_FILE.txt;
Tr = RESULT_FILE;

% extract the samples belonging to different classes
I = find(Tr(:,3) == 0);  % find the row indices for the 1st class, labeled as 0
Tr0 = Tr(I,1:2);
save Tr0;                % so that you can use it directly next time 

I = find(Tr(:,3) == 1);  % find the row indices for the 2nd class, labeled as 1
Tr1 = Tr(I,1:2);
save Tr1;

% plot the samples
plot(Tr0(:,1),Tr0(:,2),'r*'); % use "red" for class 0
hold on;           % so that the future plots can be superimposed on the previous ones
plot(Tr1(:,1),Tr1(:,2),'go'); % use "green" for class 1

%I added this part
load synth.te;
Train = synth;

I = find(Train(:,3) == 0);  % find the row indices for the 1st class, labeled as 0
Train0 = Train(I,1:2);
save Train0;                % so that you can use it directly next time 

I = find(Train(:,3) == 1);  % find the row indices for the 2nd class, labeled as 1
Train1 = Train(I,1:2);
save Train1;

% plot the samples
plot(Train0(:,1),Train0(:,2),'wo'); % use "red" for class 0
hold on;           % so that the future plots can be superimposed on the previous ones
plot(Train1(:,1),Train1(:,2),'bo'); % use "green" for class 1





