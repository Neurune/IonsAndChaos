clear all; close all; clc


P1 = importdata(['Datafiles/CheckNormSleep_Poin_cc0.txt']);
P2 = importdata(['Datafiles/CheckOnlyOne_Poin_cc0.txt']);
P3 = importdata(['Datafiles/CheckOnlySix4_Poin_cc0.txt']);
P4 = importdata(['Datafiles/CheckOnlyExt1_Poin_cc0.txt']);
P5 = importdata(['Datafiles/CheckOnlyAll_Poin_cc0.txt']);
P6 = importdata(['Datafiles/CheckSleepWithg3_Poin_cc0.txt']);
P7 = importdata(['Datafiles/CheckOnlySixWithg3_Poin_cc0.txt']);
P8 = importdata(['Datafiles/CheckOnlyThreeWithg3_Poin_cc0.txt']);
P9 = importdata(['Datafiles/CheckOnlyTenWithg3_Poin_cc0.txt']);
P10 = importdata(['Datafiles/CheckNormAwake_Poin_cc0.txt']);
N1 = length(P1(3:end,1))/100.0; N2 = mean(P1(3:end,2)); N3 = std(P1(3:end,2));
% N1 = 1.0; N2 = 1.0; N3 = 1.0;

X = [ 
      (length(P1(3:end,1))/100.0)/N1 mean(P1(3:end,2))/N2 std(P1(3:end,2))/N3;
      (length(P2(3:end,1))/100.0)/N1 mean(P2(3:end,2))/N2 std(P2(3:end,2))/N3;
      (length(P3(3:end,1))/300.0)/N1 mean(P3(3:end,2))/N2 std(P3(3:end,2))/N3;
      (length(P4(3:end,1))/100.0)/N1 mean(P4(3:end,2))/N2 std(P4(3:end,2))/N3;
      (length(P5(3:end,1))/300.0)/N1 mean(P5(3:end,2))/N2 std(P5(3:end,2))/N3;
      (length(P6(3:end,1))/100.0)/N1 mean(P6(3:end,2))/N2 std(P6(3:end,2))/N3;
      (length(P7(3:end,1))/100.0)/N1 mean(P7(3:end,2))/N2 std(P7(3:end,2))/N3;
      (length(P8(3:end,1))/100.0)/N1 mean(P8(3:end,2))/N2 std(P8(3:end,2))/N3;
      (length(P9(3:end,1))/100.0)/N1 mean(P9(3:end,2))/N2 std(P9(3:end,2))/N3;
      (length(P10(3:end,1))/100.0)/N1 mean(P10(3:end,2))/N2 std(P10(3:end,2))/N3;
    ];
c = categorical({'Sleep','I1','I6','I11','All','g = 0.9g','g = 0.9g + I6','g = 0.9g + I3','g = 0.9g + I11','Awake'});
y = bar(X); 
set(gca,'XTickLabel', {'Sleep','I1','I6','I11','All','g09','g09+I6','g09+I3','g09+I11','Awake'});
y(1).LineWidth = 2;
y(1).EdgeColor = 'black';
y(1).FaceColor = 'green';
y(2).LineWidth = 2;
y(2).EdgeColor = 'black';
y(2).FaceColor = 'red';
y(3).LineWidth = 2;
y(3).EdgeColor = 'black';
y(3).FaceColor = 'white';
goodplot; hold on

% errorbar(X,EX,'.')
