clc
clear

q = 1.602E-19;
mo = 9.11e-31;
mn = 0.26*mo;
kb = 1.381e-23;
T = 300;
Pscat = 1-exp(-1e-14/0.2e-12);
Vx = 0.5;
Vy = 0.2;
Efieldx = Vx/2E-7; %E=v/m
Forcex = Efieldx*q;
Efieldy = Vy/1E-7; %E=v/m
Forcey = Efieldy*q;
Vavg = 0;


%Acceleration Due to Efield
accelx = Forcex/mo;
accely = Forcey/mo;

%Initialise the particles
initialX = 200e-9*rand(10,1);
initialY = 100e-9*rand(10,1);
axis ([0 200e-9 0 100e-9])

%Initialise angles
angleRad = 2*pi*rand(10,1);

%Set velocity
vth = sqrt((kb*T)/mn);

%Maxwell Boltzman Inital Velocity
MD1 = randn(10,1).*(vth/sqrt(2));
MD2 = randn(10,1).*(vth/sqrt(2));
MaxwellBoltzman = sqrt((MD1).^2+(MD2).^2);
initialRV = MaxwellBoltzman;

% figure(1)
% velocity = histogram(initialRV,20);

scat = 0;


velocityX = initialRV.*cos(angleRad);
velocityY = initialRV.*sin(angleRad);

newX = initialX + velocityX*1e-14;
newY = initialY + velocityY*1e-14;


for time = 0:1e-14:3e-13
   
    
    %Check for Scatter
    Escat = rand(10,1) < Pscat;
     if newX(Escat) > 0
        %Rethermalize
        MD1 = randn(10,1).*(vth/sqrt(2));
        MD2 = randn(10,1).*(vth/sqrt(2));
        MaxwellBoltzman = sqrt((MD1).^2+(MD2).^2);
        initialRV = MaxwellBoltzman;

        %Find New Velocities
%         voldX = velocityx;
%         vnewX = (accel)*1E-14 + voldX;
        angleRad = 2*pi*rand(10,1);
        velocityX(Escat) = initialRV(Escat).*cos(angleRad(Escat)) + (accelx*1E-14);
        velocityY(Escat) = initialRV(Escat).*sin(angleRad(Escat)) + (accely*1E-14);
        
         %Find new positions
         newX = initialX + velocityX*1e-14;
         newY = initialY + velocityY*1e-14;
        
%         figure(1)
%         velocity = histogram(initialRV,25);

        %Mean Free Path/Time Between Collisions
        Vavgold = Vavg;
        scat = scat+ sum(Escat);
        tauMN = (time*10)/scat;
        Vavg = mean((velocityX.^2) + (velocityY.^2));
       
        MFP = tauMN*Vavg;
        

%  
%         figure(2)
%         subplot(2,1,1)
%         title('Mean Free Time')
%         xlabel('Total Time in Sim (s)')
%         ylabel('Mean Free Time (s)')
%         plot(time, tauMN, 'b.')
%         hold on
%         
%         figure(2)
%         subplot(2,1,2)
%         xlabel('Total Time in Sim (s)')
%         ylabel('Mean Free Path (m)')
%         title('Mean Free Path')
%         plot(time, MFP, 'g.')
%         hold on
        
     else
        Vavgold = Vavg;
        Vavg = mean((velocityX.^2) + (velocityY.^2));
      
        velocityX = velocityX + accelx*1e-14;
        velocityY = velocityY + accely*1e-14;
         %Find new positions
         newX = initialX + velocityX*1e-14 + (accelx*(1E-14^2));
         newY = initialY + velocityY*1e-14 + (accely*(1E-14^2));
        
         
     end
     
%      Vavgold = Vavg;
     
    %Track Average Velocity
    figure(1)
    title('Average Velocity')
    xlabel('Total Time in Sim (s)')
    ylabel('Average Particle Velicity (m/s)')
    plot(time, Vavg, 'b.')
%     plot([time-1e-14 time], [Vavgold Vavg], 'b')
    hold on
     
    
    %Find Drift Current Density
    PDCX = velocityX >= 0;
    driftElecX = sum(velocityX(PDCX))*q;
    Jx = driftElecX/(1e-7*2e-7); 
%     figure(1)
%     title('Current Density in X')
%     plot(time, Jx)
    
    PDCY = velocityY >= 0;
    driftElecY = sum(velocityY(PDCY))*q;
    Jy = driftElecY/(2e-7*1e-7); 
%     figure(2)
%     title('Current Density in Y')
%     plot(time, Jy)
    
    
    figure(2)
    plot(time, Jx, 'r.', time, Jy, 'b.') 
    hold on
    
    %Find temperature
    Vavg = mean((velocityX.^2) + (velocityY.^2));
    T = (mn*Vavg)/(kb);

    
    %Check X boundary conditions
    [NH,IH] = max(newX);
    [NL,IL] = min(newX);
    
    upperX = newX >= 200e-9;
    newX(upperX)= newX(upperX)-200e-9;
    initialX(upperX) = newX(upperX);
    
    lowX = newX <= 0;
    newX(lowX) = newX(lowX)+200e-9;
    initialX(lowX) = newX(lowX);
        
    
    
    %Check Y boundary conditions
    [NumH,IndexH] = max(newY);
    [NumL,IndexL] = min(newY);
    
    upperY = newY >= 100e-9;
    velocityY(upperY)= -velocityY(upperY);
    
    lowY = newY <= 0;
    velocityY(lowY) = -velocityY(lowY);
    
    
   
    
    xp = [initialX newX];
    yp = [initialY newY];
    c =  lines(10);

    figure(3)
    title('Confined Atoms');
    P = plot(xp.', yp.');
    set(P, {'color'}, num2cell(c,2));    
    hold on

    axis ([0 200e-9 0 100e-9])
    
    
    figure(4)
    title('Average Temperature (C)');
    plot(time, T, 'r.')
    hold on  
    
    
    
    
    
    
%     %Density Plots
% 
%     %Find Temperature per Atom
%     Vatoms = velocityX.^2 + velocityY.^2;
%     tempAtoms = (mn.*Vatoms)./kb;
%     
%     %Make Plot
%     x = linspace(0,2e-7,10);
%     y = linspace(0,1e-7,10);
%     [Xt,Yt] = meshgrid(newX,newY);
%     Z = tempAtoms;
%     surf(Xt,Yt,Z)
%        
%         
%         
% %     pos = meshgrid(newX, newY);
% %     map = meshgrid(pos, substrate);        
% %     
% %     x = linspace(0, 2e-7,100);
% %     y = linspace(0, 1e-7,100);
% %     [X,Y] = meshgrid(x,y);
% %     
% %     for i = 0:x
% %         for j = 0:y
% %             if particle
%     
%     figure(5)
%     title('Electron Temperature')
%     surf(tempPlot)

 




    initialX = newX;
    initialY = newY;
    
    
    pause(0.001)
    time

end

del = 0.05e-7;
nx = 40;
ny = 20;
ED = zeros(ny,nx);
Temp = zeros(ny,nx);

for i = 1:nx
    for j = 1:ny
        LowerBoxLeftX1 = newX < (i*del);
        LowerBoxLeftX2 = newX > (i*del - del); 
        LowerBoxLeftX = LowerBoxLeftX1>0 & LowerBoxLeftX2>0;      
        
        
        LowerBoxLeftY1 = newY < (j*del);
        LowerBoxLeftY2 = newY > (j*del - del);
        LowerBoxLeftY = LowerBoxLeftY1>0 & LowerBoxLeftY2>0;
        
        particles = LowerBoxLeftX>0 & LowerBoxLeftY>0;
        
        Xvelocity = velocityX(particles);
        Yvelocity = velocityY(particles);

        
        TotalElectrons = LowerBoxLeftX>0 & LowerBoxLeftY>0;
        Vavg = mean((Xvelocity.^2) + (Yvelocity.^2));
        
        
        Temp(j,i) = (mn*Vavg)/(kb);
        ED(j,i) = sum(TotalElectrons);
       
    end
end

x = linspace(0,2e-7,40);
y = linspace(0,1e-7,20);
[X,Y] = meshgrid(x,y);
[Xp,Yp] = meshgrid(newX,newY);

figure(5)
title('Electron Density')
xlabel('X Distance In Plane (1square = 0.05e-7m')
ylabel('Y Distance In Plane (1square = 0.05e-7m')
surf(X,Y,ED)


notTemp = isnan(Temp);
Temp(notTemp) = 0;

figure(7)
title('Temperature Map')
xlabel('X Distance In Plane (1square = 0.05e-7m')
ylabel('Y Distance In Plane (1square = 0.05e-7m')
surf(X,Y,Temp)







hold off

