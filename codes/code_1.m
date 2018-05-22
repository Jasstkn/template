% ns_Izh_006.m
%24 nov 2016


close all
clear all
clc

% Стандартные параметры --------------------------------------------------------
%  C = 100; vr = -60; vt = -40; k = 0.7; % параметры для нейрона типа RS
%  a = 0.03; b = -2; c = -50; d = 100;
%  vPeak = 35;        % сброс спайка  [мВ]
%  dt = 1;            % временной шаг  [мс]

% Входные данные ---------------------------------------------------------------
flagI = 1;           % моделирование требуемого нейрона
flagS = 4;           % flags = 5 for spike rates
Sstep = 200;          % шаг для входного сигнала
S1 = 0.1; S2 = 0.9;  % проценты для тока включения и выключения
dt = 0.5;            
NT = 800;              % количество временных итераций

switch flagI
   case 1   % default
       C = 100; vr = -60; vt = -40; k = 0.7; 
       a = 0.03; b = -2; c = -50; d = 100; 
       vPeak = 35;        

   case 2   % miscellaneous
       C = 50; vr = -60; vt = -40; k = 0.7; 
       a = 0.01; b = 5; c = -55; d = 150; 
       vPeak = 10;        
    
    case 3   %   Intrinsically Bursting (IB) Neurons  
        C = 150; vr = -75; vt = -45; k = 1.2; 
        a = 0.01; b = 5; c = -56; d = 130; 
        vPeak = 50; 
        
    case 4   %   Chattering (CH) neuronsdefault
       C = 50; vr = -60; vt = -40; k = 1.5; 
       a = 0.03; b = 1; c = -40; d = 150; 
       vPeak = 35;        
    


end  % switch


% Установка требуемых параметров
  t = 0:dt:(NT-1) * dt;   % время [мс]

   v = vr*ones(1,NT); u = 0*v; % инициирование переменной 
   S = zeros(1,NT);            % стимул

   S(round(S1*NT):round(S2*NT)) =  Sstep ;   % шаговая функция [70]
   %S = 0.411 .* t;                         % ramp  flagS = 5          
   %S(round(0.25*S1*NT):round(0.75*S2*NT)) = 70;


% Решение прямым методом Эйлера
for m = 1:NT-1  
  vT     = v(m)+ (dt/2) * (k*(v(m) - vr)*(v(m) - vt)-u(m) + S(m))/C;
  v(m+1) = vT + (dt/2)  * (k*(v(m) - vr)*(v(m) - vt)-u(m) + S(m))/C;
 % u(m+1) = u(m) + dt * a*(b*(v(m)-vr)-u(m));
   u(m+1) = u(m) + dt * a*(b*(v(m+1)-vr)-u(m));
    if v(m+1)>= vPeak  % условие появления спайка
       v(m)   = vPeak; % Амплитуда спайкового сигнала
       v(m+1) = c; % сброс мембранного потенциала
       u(m+1) = u(m+1) + d; % обновление переменной восстановления
   end; % if
   
end;

out1 = ones(1,NT);
for m = 1:NT-1
    if v(m) >= 35
        out1(m) = 1;
    else
        out1(m) = 0;
    end;
end;


% Firing rates: flagF = 1 (yes) flagF = 0 (no) ===========================
   flagF = 1; 
   if flagF == 1
      indexFire = find(v == vPeak);           % индексы для спайков
      indexFire2 = indexFire(2:end);
      indexFire1 = indexFire(1:end-1);
      ISI = dt .* (indexFire2 - indexFire1);  % Межспайковый интервал
      fireRate = 1000 ./ISI;                  % скорость возникновения спайкового сигнала
      f = mean(fireRate);                     % средняя скорость
      
      F = zeros(1,length(indexFire1));        % average fire rate  [Hz]
      for cc = 2:length(indexFire1)
         F(cc) = (fireRate(cc)+fireRate(cc-1))/2;
      end  % for
      
   disp('Межспайковый период   [мс]  ');
   fprintf(' %2.2f ',ISI);
   disp('   ');
   disp('Скорость возникновения спайкового сигнала  [Гц]   ');
   fprintf(' %2.2f ',fireRate);
   disp('  ');
   fprintf('Средняя скорость  f =  %2.2f \n',f);
   disp('  ');
   
   if flagS == 5     % ramp input stimuls
      figure(5)
      fs = 12;
      set(gcf,'units','normalized');
      set(gcf,'position',[0.42 0.05 0.2 0.2]);
      xP = S(indexFire1);  yP = F;
      plot(xP,yP,'b','linewidth',2);
      ylabel('f  [Гц]','fontsize',fs);
      xlabel('S  [пА]','fontsize',fs);
      axis([0 max(xP) 0 1.2*max(yP)])
      grid on
      set(gca,'fontsize',fs);
   end  % if
   
   end  % flagF


% Построение графиков   
subplot(2,2,1);
fs = 12;
set(gcf,'units','normalized');
plot(t, v,'lineWidth',2); 
xlabel('t [мс]','fontsize',fs);
ylabel('Мембранный потенциал [мВ]','fontsize',fs);
grid on
set(gca,'fontsize',fs);
axis([0 1.01*max(t) 1.1*min(v) 1.2*max(v)])
title('Стимул')

subplot(2,2,2);
   fs = 12;
   set(gcf,'units','normalized', 'color','w');
   plot(t, S,'lineWidth',2); 
   xlabel('t [мс]','fontsize',fs);
   ylabel('I [пA]','fontsize',fs);
   ylim([0 220])
   grid on
   set(gca,'fontsize',fs);
   title('Мембранный потенциал')

subplot(2,2,[3,4]);
   fs = 12;
   set(gcf,'units','normalized');
   plot(t, out1,'lineWidth',2); 
   xlabel('t [мс]','fontsize',fs);
   ylabel('Выходной сигнал','fontsize',fs);
   ylim([0 2])
   grid on
   set(gca,'fontsize',fs);
   title('Выходной сигнал')
   
suptitle('Модель Ижикевича')
   