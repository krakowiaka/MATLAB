%% Zadanie 1
clear;
x=linspace(0,1,1000);
%%Tx
figure('Name','T');
T= - 4.*x.^4.* tan(x.^4) -x;
plot(x,T);
hold on;
xlabel('x');
ylabel('T(x)');
title('Wykres 1.: Wspólczynnik T(x)');
grid on;

%%K1
figure('Name','K1');
K1= - x.^4.* tan(x.^4);
plot(x,K1);
hold on;
xlabel('x');
ylabel('K1(x)');
title('Wykres 2.: Wspólczynnik K1(x)');
grid on;
 
%%K2
figure('Name','K2');
K2=1*ones(length(x));
plot(x,K2);
hold on;
xlabel('x');
ylabel('K2(x)');
title('Wykres 3.: Wspólczynnik K2(x)');
grid on;

%%K3
figure('Name','K3');
K3= -x-2;
plot(x,K3);
hold on;
xlabel('x');
ylabel('K3(x)');
title('Wykres 4.: Wspólczynnik K3(x)');
grid on;
 
%%K4
figure('Name','K4');
K4= -1*ones(length(x));
plot(x,K4);
hold on;
xlabel('x');
ylabel('K4(x)');
title('Wykres 5.: Wspólczynnik K4(x)');
grid on;
 
%%K5
figure('Name','K5');
K5= 1*ones(length(x));
plot(x,K5);
hold on;
xlabel('x');
ylabel('K5(x)');
title('Wykres 6.: Wspólczynnik K5(x)');
grid on;
 
 
%% Zadanie 2
bl_cal_met_maks_sum_mod=max(max((abs(T)+abs(K1)+abs(K2)+abs(K3)+abs(K4)+abs(K5)))*5*10^(-13));
bl_zad_2 = vpa(bl_cal_met_maks_sum_mod);

%%Zadanie 3
syms x_3;
syms blad [6 1];
f = (cos(((x_3.*(1+blad(1)))^(4))*(1+blad(2))).*(1+blad(3))./((exp((x_3.*(1+blad(1))+2).*(1+blad(4)))).*(1+blad(5)))).*(1+blad(6));
blad_eps = de2bi(0:63);   
blad_eps(blad_eps==0)=-1;  %zamiana w kodzie binarnym 0 na -1
bledziki = blad_eps * 5*10^(-13);
y_3 = matlabFunction(f);

max_blad_3=0;

for n=1:64 
    
    b_wzgledny = ((y_3(bledziki(n,1),bledziki(n,2),bledziki(n,3),bledziki(n,4),bledziki(n,5),bledziki(n,6),x))-y_3(0,0,0,0,0,0,x))./(y_3(0,0,0,0,0,0,x));
    bl_koncowy = max(max(abs(b_wzgledny)));
    
    if(bl_koncowy>max_blad_3)
        max_blad_3=bl_koncowy;
    end
end

bl_met_sym = max_blad_3;
bl_zad_3 = vpa(max_blad_3);

%Zadanie 4
A = [ 58 47 -150 12 153; 144 131 -405 36 414 ; 90 80 -250 20 260 ;16 14 -45 9 46; 26 24 -75 4 81];
xT = [1,0,1,-1,1];
x_do_b = xT';
b = A*x_do_b;

epsik = 5*10^(-13);

max_blad_4 = 0;

for n=0:(2^(25)-1)
    
    blad_epsik = de2bi(n,25);
    macierz_epsik = reshape(blad_epsik,[5,5]);
    macierz_epsik(macierz_epsik==0)= -1;
    macierz_eps=macierz_epsik * epsik;
    
    A_zab = A.*(1+macierz_eps);
    x_blad = A_zab\b;
    
    bl_wzgledny_4 = norm((x_blad - x_do_b))/norm(x_do_b);
  
    if(bl_wzgledny_4>max_blad_4)
        max_blad_4=bl_wzgledny_4;
    end
end

blad_zad_4 =vpa(max_blad_4);

uwa_wsk = cond(A);

