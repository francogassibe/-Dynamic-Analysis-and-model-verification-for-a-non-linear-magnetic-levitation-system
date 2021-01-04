close all
clc,clear

% DISCRETIZACION TIEMPO Y FRECUENCIA

fs=10000;
Dt=1/fs; 
Ttol=10; 
n=1+(Ttol-Dt)/Dt;
t=0:Dt:(Ttol-Dt);


%PARAMETROS SISTEMA

m=0.09;
c=17000000;
pext=0;
P=m*9.8;

k=0.065^2*P;%el iman grande

p=P+pext;
xc=(k/p)^(1/2);


% CONDICIONES INICIALES

%x0=xc*1.01;%comportamiento lineal
x0=xc*4 %comportamiento no lineal

xd0=0;
xdd0=((p/m) + xd0*(c/m) + (k/m)*((x0))); %DESPEJADO DE EDO, PARA APLICAR EN DIF CENTRAL



% INICIALIZACION DE DIF CENTRAL

xant = x0 + Dt*xd0 + 0.5*(Dt^2)*xdd0; %Taylor
a=Ttol/Dt;
x=zeros(1,a);
x(1,1)=xant;
xmin=0.0001;
xk=x0;



% DIFERENCIA CENTRAL 
for j=2:(Ttol/Dt)

% FORMULA DE RECURRENCIA
x(1,j)=(-p + (k/((xk)^2) + (2*m*xk/(Dt^2)) - ((m/(Dt^2))-(c/2*Dt))*xant))/((c/2*Dt)+(m/(Dt^2)));

if(x(1,j)< xmin)
    x(1,j)=0;
    break;    
end

% ACTUALIZACION DE VARIABLES
xant=xk;
xk=x(1,j);
end

figure(1)
hold on;
grid on;
plot(t,x(1,:),'r');
ylabel('Posición [m]'); 
xlabel('Tiempo [t]');

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
 plot(t,acel,'y');
 

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
 title('Espectro en frecuencia');
xlabel('Frecuencia'); 
ylabel('Amplitud');
 maximo=max(abs(y))




 
 