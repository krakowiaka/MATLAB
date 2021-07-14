%%Zadanie 1

%macierz wspołczynników wyjsciowej funkcji w(x)
A = [-40180, -6272, 1751, 226, -68, -2, 1];

%macierz wspolczynników do obliczenia idealnych rozwiazan  
A_1 = [1, -2, -68, 226, 1751, -6272, -40180];

%rysowanie funkcji w(x)
x = linspace(-10,10,10000);
figure(1);
plot(x,polyval(A_1,x));
xlabel('x');
ylabel('w(x)');
title('Wykres funkcji w(x)')

%miecierz z idealnymi rozwiązaniami
wektor_rozwiazan = roots(A_1);

dokladnosc = 10^(-4);
xi_MI = -3.5 +1.5i;
xi_sie = 6.5;
ximinus1_sie = 8;
xi_MII = 4.5 + 3.5i;
ximinus1_MII = 5 + 5.5i;
ximinus2_MII = 5 + 6i;
xi_st = -6.5;

pier_zespol1 = MullerI(A_1, xi_MI, dokladnosc);
pier_zespol1_sprzezone = conj(pier_zespol1);

A_2 = def_kwa(A, pier_zespol1);
A_2_od = fliplr(A_2);

pier_rzeczy1 = sieczne(A_2_od, xi_sie, ximinus1_sie, dokladnosc);

A_3 = def_lini(A_2, pier_rzeczy1);
A_3_od = fliplr(A_3);

pier_zespol2 = MullerII(A_3_od, xi_MII, ximinus1_MII, ximinus2_MII, dokladnosc);
pier_zespol2_sprzezone = conj(pier_zespol2);

A_4 = def_kwa(A_3, pier_zespol2);
A_4_od = fliplr(A_4);

pier_rzeczy2 = styczne(A_4_od, xi_st, dokladnosc);

%%Zadanie 2

delta_w = 10.^(-4:-1:-16);

Wyniki = [1;2];
Wyniki_nie = [1;2];

for i=1:length(delta_w)
    
    pier_zespol1_2 = MullerI(A_1, xi_MI, delta_w(i));
    pier_zespol1_2_sprzezone = conj(pier_zespol1_2);

    A_2_2 = def_kwa(A, pier_zespol1_2);
    A_2_od_2 = fliplr(A_2_2);

    pier_rzeczy1_2 = sieczne(A_2_od_2, xi_sie, ximinus1_sie, delta_w(i));

    A_3_2 = def_lini(A_2_2, pier_rzeczy1_2);
    A_3_od_2 = fliplr(A_3_2);

    pier_zespol2_2 = MullerII(A_3_od_2, xi_MII, ximinus1_MII, ximinus2_MII, delta_w(i));
    pier_zespol2_2_sprzezone = conj(pier_zespol2_2);

    A_4_2 = def_kwa(A_3_2, pier_zespol2_2);
    A_4_od_2 = fliplr(A_4_2);

    pier_rzeczy2_2 = styczne(A_4_od_2, xi_st, delta_w(i));
    
    wektor_zad2 = [ pier_rzeczy1_2; pier_zespol2_2; pier_zespol2_2_sprzezone; pier_rzeczy2_2; pier_zespol1_2; pier_zespol1_2_sprzezone];

    blad_wzgledny = norm(wektor_zad2 - wektor_rozwiazan)/norm(wektor_rozwiazan);
    blad_wzgledny_nie = norm(wektor_zad2 - wektor_rozwiazan,Inf)/norm(wektor_rozwiazan,Inf);
    
    if( isnan(blad_wzgledny) || isnan(blad_wzgledny_nie))
        blad_wzgledny = eps;
        blad_wzgledny_nie = eps;
    end

    Wyniki = [Wyniki(1,:) delta_w(i);
    Wyniki(2,:) blad_wzgledny];

    Wyniki_nie = [Wyniki_nie(1,:) delta_w(i);
    Wyniki_nie(2,:) blad_wzgledny_nie];
    
end

delta_w_2 =10.^(-4:-1:-16);
figure(2)
loglog(delta_w_2,Wyniki(2,2:end),'k*');
hold on
title('Zalezność bledów względnych wektor x od parametru ∆w - norma euklidesowa')
xlabel('∆w')
ylabel('zagregowany blad')
axis([10^-16 10^-4 10^-16 10^-1])
hold off

figure(3)
loglog(delta_w_2,Wyniki_nie(2,2:end),'k*')
hold on
title('Zalezność bledów względnych wektor x od parametru ∆w - norma nieskończoności')
xlabel('∆w')
ylabel('zagregowany blad')
axis([10^-16 10^-4 10^-16 10^-1])
hold off
 
function y = MullerI(w, xi_MI, dokladnosc)

w1 = polyder(w);
w2 = polyder(w1);

iter_MI = 0;
while iter_MI<1000
    a_MI = 0.5 * polyval(w2, xi_MI);
    b_MI = polyval(w1, xi_MI);
    c_MI = polyval(w, xi_MI);
    xiplus1_MI = xi_MI - (2*c_MI)/(b_MI +sign(b_MI) * sqrt(b_MI^2 - 4*a_MI*c_MI));
    xi_MI = xiplus1_MI;
    iter_MI = iter_MI+1;
    if (abs(polyval(w,xi_MI)) < dokladnosc)
       break;
    end
end

y = xi_MI;
end

function y = sieczne(w, xi_sie,ximinus1_sie, dokladnosc)

iter_sie =0;
while iter_sie < 1000
    iter_sie = iter_sie +1;
    xiplus1_sie = xi_sie - ((xi_sie - ximinus1_sie)/(polyval(w,xi_sie)-polyval(w,ximinus1_sie)))*polyval(w,xi_sie);
    ximinus1_sie = xi_sie;
    xi_sie = xiplus1_sie;
    if(abs(polyval(w,xi_sie)) < dokladnosc)
        break;
    end
end

y = xi_sie;
end

function y = MullerII(w, xi_MII, ximinus1_MII, ximinus2_MII, dokladnosc)

iter_MII=0;
while iter_MII < 1000
    macierz_1 = [ -(xi_MII - ximinus1_MII)^2, (xi_MII - ximinus1_MII); 
        -(xi_MII - ximinus2_MII)^2,(xi_MII - ximinus2_MII)];
    macierz_2 = [polyval(w,xi_MII) - polyval(w,ximinus1_MII);
        polyval(w,xi_MII) - polyval(w,ximinus2_MII)];
    w_ai_bi = macierz_1\macierz_2;
    ai_II = w_ai_bi(1,1);
    bi_II = w_ai_bi(2,1);
    ci_II = polyval(w,xi_MII);
    xiplus1_MII = xi_MII - (2*ci_II)/(bi_II + sign(bi_II)*sqrt(bi_II-4*ai_II*ci_II));
    ximinus2_MII = ximinus1_MII;
    ximinus1_MII = xi_MII;
    xi_MII = xiplus1_MII;
    iter_MII =iter_MII + 1;
    
    if (abs(polyval(w,xi_MII) == 0 || isnan(polyval(w,xi_MII))))
       break;
    end
    if(abs(polyval(w,xi_MII)) < dokladnosc)
        break;
    end
end

y = xi_MII;
end

function y = styczne(w, xi_st, dokladnosc)

w1=polyder(w);

iter_st = 0;
while iter_st < 1000
    iter_st = iter_st+1;
    xiplus1_st = xi_st - (polyval(w,xi_st)/polyval(w1,xi_st));
    xi_st = xiplus1_st;
    if(abs(polyval(w,xi_st)) < dokladnosc)
        break;
    end
end

y = xi_st;
end

function y = def_kwa(A, pierw_zesp)

N = size(A,2);
p = 2*real(pierw_zesp);
r = -(abs(pierw_zesp)^2);
B(N-2) = A(N);
B(N-3) = A(N-1)+p*B(N-2);
for n=N-4:-1:1
    B(n) = A(n+2) +p*B(n+1) + r*B(n+2);
end

y = B;
end

function y = def_lini(A, pierw_rzecz)

N = size(A,2);
B(N-1) = A(N);
for n=N-2:-1:1
    B(n) = A(n+1)+pierw_rzecz.*B(n+1);
end

y = B;
end
