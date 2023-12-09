%Filtrado lineal en el dominio frecuencial
clear all
clc

N=128;
t= 0:0.01:N-1;
xt=2*cos(2*pi*1/32*t)-3*sin(2*pi*1/16*t)+4*cos(2*pi*1/8*t)+1*cos(2*pi*1/16*t);
figure('Name','señal x(t)')
plot(t,xt)
n= 0:N-1;
xn=2*cos(2*pi*1/32*n)-3*sin(2*pi*1/16*n)+4*cos(2*pi*1/8*n)+1*cos(2*pi*1/16*n);
xn2=0*cos(2*pi*1/32*n)-3*sin(2*pi*1/16*n)+0*cos(2*pi*1/8*n)+1*cos(2*pi*1/16*n);
figure('Name','señal x(n)')
stem(n,xn)

%Espectro de x(n) usando fft
xk=fft(xn);
%Filtrado pasabajos ideal wc=1/10
wc=1/10;
wct=floor(wc*N);
%Filtrado pasaltos ideal wc=1/20
wc1=1/20;
wct1=floor(wc*N);
xkf=xk;
for i=1:N
    if (i>=wct) && (i<N-wct)
        xkf(i)=0;
    elseif (i<=wct1) && (i<N-wct1)
        xkf(i)=0;
    end
end
figure('Name','Espectro (Re)')
stem(n,real(xk),'bd','LineWidth',1.5)
hold on
stem(n,real(xkf),'r*')
grid on
legend('original','filtrada')

figure('Name','Espectro (Im)')
stem(n,imag(xk),'bd','LineWidth',1.5)
hold on
stem(n,imag(xkf),'r*')
grid on
legend('original','filtrada')

figure('Name','Espectro (magnitud)')
stem(n,abs(xk),'bd','LineWidth',1.5)
hold on
stem(n,abs(xkf),'r*')
grid on
legend('original','filtrada')

figure('Name','Espectro (Fase)')
stem(n,angle(xk),'bd','LineWidth',1.5)
hold on
stem(n,angle(xkf),'r*')
grid on
legend('original','filtrada')


% Reconstrucción señal filtrada
xnf=real(ifft(xkf));
figure('Name','Señal x(n) ');
stem(n,xn,'bd','LineWidth',1.5)
hold on
stem(n,xnf,'r*')
stem(n,xn2,'y*')
legend('original','filtrada','1/16')

xtf = ifft(xnf);
figure('Name','Señal x(t) filtrada ');
plot(t,xt)
hold on
plot(xnf)
legend('original','filtrada')

