load('dane.mat');

A = mammographicmasses;
TF = ismissing(A);
T3 = rmmissing(A);
A = table2array(A);

column6 = A(:,6);

%%złośliwe
col1a = A(:,1);
col1a(column6==0) = [];
col1a(col1a == 0) = NaN;
col1a(col1a > 5) = NaN;
c1a = fillmissing(col1a, 'constant', median(col1a, 'omitnan'));

col2a = A(:,2);
col2a(column6==0) = [];
c2a = fillmissing(col2a, 'constant', mean(col2a, 'omitnan'));

col3a = A(:,3);
col3a(column6==0) = [];
c3a = fillmissing(col3a, 'constant', median(col3a, 'omitnan'));

col4a = A(:,4);
col4a(column6==0) = [];
c4a = fillmissing(col4a, 'constant', median(col4a, 'omitnan'));

col5a = A(:,5);
col5a(column6==0) = [];
c5a = fillmissing(col5a, 'constant', median(col5a, 'omitnan'));

col6a = A(:, 6);
col6a(column6==0) = [];

%%łagodne
col1b = A(:,1);
col1b(column6==1) = [];
col1b(col1b == 0) = NaN;
col1b(col1b > 5) = NaN;
c1b = fillmissing(col1b, 'constant', median(col1b, 'omitnan'));

col2b = A(:,2);
col2b(column6==1) = [];
c2b = fillmissing(col2b, 'constant', mean(col2b, 'omitnan'));

col3b = A(:,3);
col3b(column6==1) = [];
c3b = fillmissing(col3b, 'constant', median(col3b, 'omitnan'));

col4b = A(:,4);
col4b(column6==1) = [];
c4b = fillmissing(col4b, 'constant', median(col4b, 'omitnan'));

col5b = A(:,5);
col5b(column6==1) = [];
c5b = fillmissing(col5b, 'constant', median(col5b, 'omitnan'));

col6b = A(:, 6);
col6b(column6==1) = [];

c11 = cat(1,c1a,c1b);
c22 = cat(1,c2a,c2b);
c33 = cat(1,c3a,c3b);
c44 = cat(1,c4a,c4b);
c55 = cat(1,c5a,c5b);
c66 = cat(1,col6a,col6b);

N = A;
N(:,1) = c11;
N(:,2) = c22;
N(:,3) = c33;
N(:,4) = c44;
N(:,5) = c55;
N(:,6) = c66;

%standaryzacja
srednia = mean(N);
odchylenie = std(N);

pomoc1 = size(N);
daneStand = zeros(pomoc1);

for i = 1:pomoc1
    daneStand(i,1) = (N(i,1)-srednia(1,1))./odchylenie(1,1);
    daneStand(i,2) = (N(i,2)-srednia(1,2))./odchylenie(1,2);
    daneStand(i,3) = (N(i,3)-srednia(1,3))./odchylenie(1,3);
    daneStand(i,4) = (N(i,4)-srednia(1,4))./odchylenie(1,4);
    daneStand(i,5) = (N(i,5)-srednia(1,5))./odchylenie(1,5);
    daneStand(i,6) = N(i, 6);
end


%podział na dane testowe i treningowe
cv = cvpartition(size(daneStand,1),'HoldOut',0.15);
idx = cv.test;
dataTrain = daneStand(~idx,:);
dataTest  = daneStand(idx,:);

dataTrain_test = dataTrain; %pełen zbiór treningowy
dataTrain_ocz = (dataTrain(:,6))'; %wartosci oczekiwane zbioru treningowego (wyjscie)
dataTrain(:,6) = [];  %zbiór treningowy bez oczekiwanych wartości (wejscie)

dataTest_test = dataTest; %pełen zbiór testowy
dataTest_ocz = (dataTest(:,6))'; %wartosci oczekiwane zbioru testowego (wyjscie)
dataTest(:,6) = []; %zbiór testowy bez oczekiwanych wartości (wejscie)


W_ww = 2*rand(3,5)-1;
W_wu = 2*rand(1,3)-1;
C = 100000;
blad_mse_train = 0;
Err = zeros(C,3);
Tre = size(dataTrain, 1);

for cykl = 1:C
      [W_ww, W_wu, E] = backProp(W_ww, W_wu, dataTrain, dataTrain_ocz);
      Err(cykl, 1) = cykl;
      Err(cykl, 2) = E;
      blad_mse_train = E*2/Tre;
      Err(cykl, 3) = blad_mse_train;
end

figure();
semilogx(Err(:,1),Err(:,2));
title('Wykres wartosci błędu w funkcji iteracji train');
xlabel('numer iteracji');
ylabel('wartość błędu');

Test = size(dataTest_test,1);
blad_test = 0;
  for q = 1:Test
      dataTest1 = dataTest(q,:)';
      x_test = W_ww*dataTest1;
      y_test = tangh(x_test);
      x_test1 = W_wu*y_test;
      y_test1 = tangh(x_test1);
      
      if y_test1 >= 0.5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

            y_test1 = 1;
      else
          y_test1 = 0;
      end
      
      Y_test(q,1) = y_test1;
      e_test = dataTest_ocz(1,q)-y_test1;
      blad_test = blad_test+0.5*e_test^2;
  end

blad_mse_test1 = blad_test*2/Test;


TN = 0;
TP = 0;
FP = 0;
FN = 0;

for m=1:Test
   
    if Y_test(m,1) == 0 && dataTest_ocz(m) == 0 
        TN = TN + 1;
    elseif Y_test(m,1) == 1 && dataTest_ocz(m) == 1
        TP = TP + 1;
    elseif Y_test(m,1) == 0 && dataTest_ocz(m) == 1
        FN = FN + 1;
    else
        FP = FP + 1;
    end
end

%czułość

TPR = TP/(TP+FN);

%specyficzność

TNR = TN/(FP+TN);

%% 


%implementacja funkcji aktywacji wartstwy wyjściowej i ukrytej
%funkcja tangens hiperboliczny
function [y] = tangh(v)
    y = tanh(v);
end

%pochodna
 function [dtang] = difftanh(v)
     dtang = 1 -tanh(v).*tanh(v);
 end
 


%algorytm uczenia - propagacja wsteczna
function [W_ww, W_wu, E] = backProp(W_ww, W_wu, X, D)
    alfa = 0.001; %współczynnik szybkości uczenia
    T = size(X,1);
    E = 0;
    for i = 1:T  
        x = X(i,:)';
        d = D(i);
        v1 = W_ww*x; %obliczenie pobudzenia warstwy ukrytej
        y1 = tangh(v1); %obliczenie stanu wyjścia wartswy ukrytej
        v = W_wu*y1 ; %obliczenie pobudzenia warstwy wyjściowej
        y = tanh(v); %obliczenie stanu wyjścia warstwy wyjściowej
        
        e = d - y; %różnica między wartością oczekiwaną, a uzyskaną y
        delta = difftanh(v)*e; 
        e1 = W_wu'*delta; 
        delta1 = difftanh(v1).*e1;
        dW_ww = alfa.*delta1*x';
        W_ww = W_ww + dW_ww; 
        dW_wu = alfa*delta*y1';
        W_wu = W_wu + dW_wu; 
        E = E + 0.5*e^2;
    end
    E = E/T;
end