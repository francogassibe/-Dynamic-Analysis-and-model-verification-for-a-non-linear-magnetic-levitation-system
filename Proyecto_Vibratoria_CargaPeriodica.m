close all
clc,clear

% DISCRETIZACION TIEMPO Y FRECUENCIA

fs=10000;
Dt=1/fs; 
Ttol=13; 
n=1+(Ttol-Dt)/Dt;
t=0:Dt:(Ttol-Dt);


%PARAMETROS SISTEMA
m=0.09;
c=10000000;
P=m*9.8;
k=0.07^2*P;


frseca=3.1*10^(-3)*9.8;%%
pext=0.09*9.8*(cos(2*pi*1.4*t));
p=P+pext;
xc=(k/P)^(1/2)


% CONDICIONES INICIALES


%x0=0.2;%comportamiento no lineal
x0=xc*1.001 %comportamiento lineal
xd0=0;
xdd0=((p(1)/m) + xd0*(c/m) + (k/m)*((x0))); %DESPEJADO DE EDO, PARA APLICAR EN DIF CENTRAL



% INICIALIZACION DE DIF CENTRAL
xant = x0 + Dt*xd0 + 0.5*(Dt^2)*xdd0;
a=Ttol/Dt;
x=zeros(1,a);
x(1,1)=xant;
xmin=0.0001;
xk=x0;
dxpaso=xd0;


% DIFERENCIA CENTRAL 


for j=2:(Ttol/Dt)
    %friccion seca
if dxpaso<0
    p=P+pext-frseca;%%
elseif dxpaso>0
    p=P+pext+frseca;%%
else
    p=P+pext;
end


% FORMULA DE RECURRENCIA
x(1,j)=(-p(j) + (k/((xk)^2) + (2*m*xk/(Dt^2)) - ((m/(Dt^2))-(c/2*Dt))*xant))/((c/2*Dt)+(m/(Dt^2)));

if(x(1,j)< xmin)
    x(1,j)=0;
    break;    
end

% ACTUALIZACION DE VARIABLES
xant=xk;
xk=x(1,j);
dxpaso=(xk-xant)/Dt;
end

figure(1)
hold on;
grid on;
plot(t,x(1,:),'r');
title('Grafica 5');
xlabel('Posición'); 
ylabel('Tiempo');


%REDEFINICION
z = x(1,:);


%CALCULO DE DERIVADAS
 for p=2:length(z)-1
     ddx(p)=(1/(2*Dt))*((-z(p-1) + z(p+1)));
 end
 
 %COMPLETO VECTOR DE VELOCIDAD
 ddx(1)=xd0;
 ddx(length(z))=ddx(length(z)-1);
 
 
 figure(2)
 grid on;
 hold on;
 plot(t,ddx)
 
 
 
 %CALCULO DE ACELERACIONES
for g=2:length(z)-1
     acel(g)=(1/(2*Dt))*((-ddx(g-1) + ddx(g+1)));
end
 
 acel(1)=xdd0;
 acel(length(z))=acel(length(z)-1);

 figure(3)
 grid on;
 hold on;
 plot(t,acel);
 

%CALCULO ESTIMADO DE LA FRECUENCIAclc,clea WD

F(1)=1;
cont=0;
for u=2:length(z)
    if(((ddx(u-1)<0 && ddx(u)>0)))
        F(u)=1;
        cont=cont+1;
    else
        F(u)=0;
    end 
end

Frecuenciaprom = cont / Ttol


%Diagrama de fases
figure(6)
grid on;
plot(x(1,:),ddx);

figure(7)
grid on;
plot(ddx,acel,'g');



%  DISCRETIZACION DE FRECUENCIA PARA FOURIER
 fn=0:fs/n:fs-fs/n;
 y=fft(x(1,:));
 


 figure(4)
 grid on;
 plot(fn,abs(y)/length(t));
 maximo=max(abs(y))




 
 