
%%Zadanie 1
clear all
close all
clc

syms x
w_x = x.^6 - 2*x.^5 - 68*x.^4 + 226*x.^3 + 1751*x.^2 - 6272*x - 40180;

rl = -7;
ru = 7;

estymata = int(w_x,rl,ru);
estymata_modul = abs(estymata);
wartosc_estymaty = double(estymata_modul);

%%Zadanie 2
warunek = 10^(-6);

macierz = [1 1 0 0 0 0 0 ;
           1 4 1 0 0 0 0 ;
           1 3 3 1 0 0 0 ;
           7 32 12 32 7 0 0 ;
           19 75 50 50 75 19 0 ;
           41 216 27 272 27 216 41];
             
wynik = [1;2;3;4];

for N=1:6
    
    mianownik = sum(macierz(N,:));
    M = 1; %ilość podprzedziałów
    
    while M < 10000
        przedzial = linspace(rl,ru,M+1);
        suma_wynik_przedzial= 0; %wyliczone pole powierzchni
        for pp=1:M
            wezly = linspace(przedzial(pp), przedzial(pp+1), N+1);
            dlugosc_kroku_calkowania = (przedzial(pp+1) - przedzial(pp));
            wynik_podrzedzialu = (sum(macierz(N,1:N+1).*w(wezly)))* dlugosc_kroku_calkowania/mianownik;
            suma_wynik_przedzial = suma_wynik_przedzial + wynik_podrzedzialu;
            wyliczone_pole(N) = suma_wynik_przedzial; 
        end

        blad = abs(abs(wyliczone_pole(N)) - wartosc_estymaty);
        
        if(blad < warunek)
            break;
        end
    
        M=M+8;
    end
    
    wynik = [wynik(1,:) N;
        wynik(2,:) M;
        wynik(3,:) blad;
        wynik(4,:) dlugosc_kroku_calkowania];
end

figure(1)
semilogy(wynik(1,3:end),wynik(2,3:end),'k')
hold on
semilogy(wynik(1,3:end),wynik(2,3:end),'k*')
title('Zalezność ilości podprzedziałów od N')
xlabel('N')
ylabel('ilosc podprzedziałów')
grid on
hold off

figure(2)
semilogy(wynik(1,3:end),wynik(4,3:end),'k')
hold on
semilogy(wynik(1,3:end),wynik(4,3:end),'k*')
title('Zalezność długości kroku całkowania od N')
xlabel('N')
ylabel('dlugosc kroku calkowania')
grid on
hold off

%%Zadanie 3
N1 = 10:10:100;
N2 = 200:100:1000;
N3 = 2000:1000:10000;
N_glowne = [N1 N2 N3];

przedzial_w_x = linspace(-7,7,10000);
wartosci_w_x = w(przedzial_w_x);
min_funkcji_w_x = min(wartosci_w_x);

pole_prosto = 14*abs(min_funkcji_w_x);

wyniki_losowe = [1;2;3;4];
wyniki_row = [1;2;3;4];

%losowanie przypadkowe
for k =1:length(N_glowne)
   x_losowe = rand((N_glowne(k))^2,1)*14-7;
   y_losowe = (rand(N_glowne(k)^2,1))*(min_funkcji_w_x);
   uzysk_los = length(find(y_losowe<=(w(x_losowe))))/(N_glowne(k))^2; %ilosc pod wykresem
   uzysk_los_nad_wykresem = 1- uzysk_los; %ilosc nad wykresem
   I_wyliczone_los = uzysk_los_nad_wykresem*pole_prosto;
   blad_bezwzgledny_los = abs(I_wyliczone_los-wartosc_estymaty);
   
   wyniki_losowe = [ wyniki_losowe(1,:) N_glowne(k);
       wyniki_losowe(2,:) uzysk_los_nad_wykresem;
       wyniki_losowe(3,:) I_wyliczone_los;
       wyniki_losowe(4,:) blad_bezwzgledny_los];
end

%losowanie równomierne
for k =1:length(N_glowne)
   
   przedzial_row_x = linspace(rl,ru,N_glowne(k));
   przedzial_row_y = linspace(min_funkcji_w_x,0,N_glowne(k));
   
   [X_row, Y_row] = meshgrid(przedzial_row_x, przedzial_row_y);
   
   uzysk_row = numel(find(Y_row<=w(X_row)))/N_glowne(k)^2;
   uzysk_row_nad_wykresem = 1 - uzysk_row;  
    
   I_wyliczone_row = uzysk_row_nad_wykresem*pole_prosto;
 
   blad_bezwzgledny_row = abs(I_wyliczone_row-wartosc_estymaty);
   
   wyniki_row = [ wyniki_row(1,:) N_glowne(k);
       wyniki_row(2,:) uzysk_row_nad_wykresem;   
       wyniki_row(3,:) I_wyliczone_row;
       wyniki_row(4,:) blad_bezwzgledny_row];
end

figure(3)
semilogx(wyniki_losowe(1,2:end),wyniki_losowe(3,2:end),'k')
hold on
semilogx(wyniki_row(1,2:end),wyniki_row(3,2:end),'m')
semilogx(wyniki_losowe(1,2:end),wyniki_losowe(3,2:end),'k*')
semilogx(wyniki_row(1,2:end),wyniki_row(3,2:end),'m*')
title('Zalezność estymaty pola powierzchni od N')
xlabel('N')
ylabel('wartość estymaty')
legend ('losowanie przypadkowe', 'losowanie równomierne');
grid on
hold off

figure(4)
loglog(wyniki_losowe(1,2:end),wyniki_losowe(4,2:end),'k')
hold on
loglog(wyniki_row(1,2:end),wyniki_row(4,2:end),'m')
loglog(wyniki_losowe(1,2:end),wyniki_losowe(4,2:end),'k*')
loglog(wyniki_row(1,2:end),wyniki_row(4,2:end),'m*')
title('Zalezność błędu bezwzględnego od N')
xlabel('N')
ylabel('wartość błędu bezwzględnego')
legend ('losowanie przypadkowe', 'losowanie równomierne');
grid on
hold off

function y = w(x)
    y = x.^6 - 2*x.^5 - 68*x.^4 + 226*x.^3 + 1751*x.^2 - 6272*x - 40180;
end
