% ELEC 4700 Assignment 1
%Chantel Lepage
% Due date: Feb 3, 2019

clc
close all

%% Part1: Electron Modelling

k=1.38e-23;
T=300;
m=0.26*9.1e-31;
t=0.2e-12;


L=200e-9;
W=100e-9;
plotPop=10;
popNum = 3e4;

Vth=sqrt((k*T)/m)
tStep=W/Vth/100;
iterations = 1000;


Vx=(Vth/sqrt(2))*rand(popNum,1);
Vy=(Vth/sqrt(2))*rand(popNum,1);
Vdis=sqrt(Vx.^2+Vy.^2);
Vavg = mean(Vdis);
MFP = Vavg*t


electron= zeros(popNum,4);
traj=zeros(iterations, plotPop*2);
temp=zeros(iterations, 1);

temp(:,1)= 300;

for i=1:popNum
    theta=rand*2*pi;
    state(i,:)= [L*rand W*rand Vth*cos(theta) Vth*sin(theta)];    
end


for i=1:iterations
    state(:,1:2)=state(:,1:2)+tStep.*state(:,3:4);
    
    j=state(:,1)> L;
    state(j,1) = state(j,1)-L;
    
    j = state(:,1)< 0;
    state(j,1)=state(j,1)+L;
    
    j = state(:,2)> W;
    state(j,2)= 2*W - state(j,2);
    state(j,4)= -state(j,4);
    
    j = state(:,2) < 0;
    state(j,2) = -state(j,2);
    state(j,4) = -state(j,4);
    
    for j=1:plotPop
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    if mod(i,5)==0
        figure(1);
        subplot(2,1,1);
        hold off;
        plot(state(1:plotPop,1)./1e-9, state(1:plotPop,2)./1e-9, 'o');
        axis([0 L/1e-9 0 W/1e-9]);
        title(sprintf('Electrons with Fixed Velocity'));
        xlabel('x (nm)');
        ylabel('y (nm)');
        if i > 1
            subplot(2,1,2);
            hold off;
            plot(tStep*(0:i-1), temp(1:i));
            axis([0 tStep*iterations min(temp)*0.98 max(temp)*1.02]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        pause(0.05);
    end
end

figure(1);
subplot(2,1,1);
title(sprintf('Electrons with Fixed Velocity'));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 L/1e-9 0 W/1e-9]);
hold on;

for i=1:plotPop
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '.');
end


%% Part 2 Collisions with Mean Free Path (MFP)

pScat= 1-exp(-tStep/t)
vPDF = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/m));

for i=i:popNum
    state(i,:)= [L*rand W*rand random(vPDF) random(vPDF)];    
end

Vavg=sqrt(sum(state(:,3).^2)/popNum+sum(state(:,4).^2)/popNum);

for i=1:iterations
    state(:,1:2)=state(:,1:2)+tStep.*state(:,3:4);
    
    j=state(:,1)> L;
    state(j,1) = state(j,1)-L;
    
    j = state(:,1)< 0;
    state(j,1)=state(j,1)+L;
    
    j = state(:,2)> W;
    state(j,2)= 2*W - state(j,2);
    state(j,4)= -state(j,4);
    
    j = state(:,2) < 0;
    state(j,2) = -state(j,2);
    state(j,4) = -state(j,4);
    
    j = rand(popNum, 1) < pScat;
    state(j,3:4) = random(vPDF, [sum(j),2]);
    
    temp(i)=(sum(state(:,3).^2)+sum(state(:,4).^2))*m/k/2/popNum; 
    
    for j=1:plotPop
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    
    if mod(i,5)==0
        figure(2);
        subplot(3,1,1);
        hold off;
        plot(state(1:plotPop,1)./1e-9, state(1:plotPop,2)./1e-9, 'o');
        axis([0 L/1e-9 0 W/1e-9]);
        title(sprintf('Electrons with Collision and MFP'));
        xlabel('x (nm)');
        ylabel('y (nm)');
        if i > 1
            subplot(3,1,2);
            hold off;
            plot(tStep*(0:i-1), temp(1:i));
            axis([0 tStep*iterations min(temp)*0.98 max(temp)*1.02]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        subplot(3,1,3);
        v = sqrt(state(:,3).^2 + state(:,4).^2);
        title('Histogram of Electron Speeds');
        histogram(v);
        xlabel('Speed (m/s)');
        ylabel('Number of particles');
        pause(0.05);
    end
end

figure(2);
subplot(3,1,1);
title(sprintf('Electrons with Collision and MFP'));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 L/1e-9 0 W/1e-9]);
hold on;
for i=1:plotPop
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '.');
end

%the mean temperature of the plot is still around 300K even with the
%fluxtuations.

%% Part 3 Enchanements 

% The non-periodic top and bottom boundaries can be set to be either
% specular (1) or diffusive (0)
top_specular = 0;
bottom_specular = 0;



boxes = 1e-9.*[80 120 0 40; 80 120 60 100];
boxes_specular = [0 1];


for i = 1:popNum
    theta = rand*2*pi;
    state(i,:) = [L*rand W*rand random(vPDF) random(vPDF)];
    
    if (state(i,2)>60e-9 &(state(i,1)>80e-9 & state(i,1)<120e-9))  | (state(i,2)< 40e-9 &(state(i,1)>80e-9 & state(i,1)<120e-9))
        state(i,1:2) = [L*rand W*rand];
    end
end


for i = 1:iterations
    state(:,1:2) = state(:,1:2) + tStep.*state(:,3:4);
    
    j = state(:,1) > L;
    state(j,1) = state(j,1) - L;
    
    j = state(:,1) < 0;
    state(j,1) = state(j,1) + L;
    
    j = state(:,2) > W;

    if(top_specular)
        state(j,2) = 2*W - state(j,2);
        state(j,4) = -state(j,4);
    else 
        state(j,2) = W;
        v = sqrt(state(j,3).^2 + state(j,4).^2);
        theta = rand([sum(j),1])*2*pi;
        state(j,3) = v.*cos(theta);
        state(j,4) = -abs(v.*sin(theta));
    end
    
    j = state(:,2) < 0;
    
    if(bottom_specular)
        state(j,2) = -state(j,2);
        state(j,4) = -state(j,4);
    else 
        state(j,2) = 0;
        v = sqrt(state(j,3).^2 + state(j,4).^2);
        theta = rand([sum(j),1])*2*pi;
        state(j,3) = v.*cos(theta);
        state(j,4) = abs(v.*sin(theta));
    end

    for j=1:popNum
        if (state(j,2)>60e-9 &(state(j,1)>80e-9 & state(j,1)<120e-9)) 
            boxNum = 1;
        elseif (state(j,2)< 40e-9 &(state(j,1)>80e-9 & state(j,1)<120e-9))
                boxNum = 2;
        else 
            boxNum = 0;
        end
        while(boxNum ~= 0)
      
            xDist = 0;
            newX = 0;
            if(state(j,3) > 0)
                xDist = state(j,1) - boxes(boxNum,1);
                newX = boxes(boxNum,1);
            else
                xDist = boxes(boxNum,2) - state(j,1);
                newX = boxes(boxNum,2);
            end
            
            y_dist = 0;
            new_y = 0;
            if(state(j,4) > 0)
                y_dist = state(j,2) - boxes(boxNum, 3);
                new_y = boxes(boxNum, 3);
            else
                y_dist = boxes(boxNum, 4) - state(j,2);
                new_y = boxes(boxNum, 4);
            end
            
            if(xDist < y_dist)
                state(j,1) = newX;
                if(~boxes_specular(boxNum))
                    sgn = -sign(state(j,3));
                    v = sqrt(state(j,3).^2 + state(j,4).^2);
                    theta = rand()*2*pi;
                    state(j,3) = sgn.*abs(v.*cos(theta));
                    state(j,4) = v.*sin(theta);
                else 
                    state(j,3) = -state(j,3);
                end
            else
                state(j,2) = new_y;
                if(~boxes_specular(boxNum))
                    sgn = -sign(state(j,4));
                    v = sqrt(state(j,3).^2 + state(j,4).^2);
                    theta = rand()*2*pi;
                    state(j,3) = v.*cos(theta);
                    state(j,4) = sgn.*abs(v.*sin(theta));
                else
                    state(j,4) = -state(j,4);
                end
            end
             boxNum = 0;
        end
    end
    
    
    j = rand(popNum, 1) < pScat;
    state(j,3:4) = random(vPDF, [sum(j),2]);
    

    temp(i) = (sum(state(:,3).^2) + sum(state(:,4).^2))*m/k/2/popNum;
    

    for j=1:plotPop
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end 
    

    if  mod(i,5) == 0
        figure(3);
        subplot(3,1,1);
        hold off;
        plot(state(1:plotPop,1)./1e-9, state(1:plotPop,2)./1e-9, 'o');
        hold on;

        for j=1:size(boxes,1)
           plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)]./1e-9,[boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)]./1e-9, 'k-');
        end
        
        axis([0 L/1e-9 0 W/1e-9]);
        title(sprintf('Electrons',plotPop, popNum));
        xlabel('x (nm)');
        ylabel('y (nm)');
        if i > 1
            subplot(3,1,2);
            hold off;
            plot(tStep*(0:i-1), temp(1:i));
            axis([0 tStep*iterations min(temp(1:i))*0.98 max(temp)*1.02]);
            title('Semiconductor Temperature');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
        end
        
        subplot(3,1,3);
        v = sqrt(state(:,3).^2 + state(:,4).^2);
        title('Histogram of Electron Speeds');
        histogram(v);
        xlabel('Speed (m/s)');
        ylabel('Number of particles');
        
        pause(0.05);
    end
end


figure(3);
subplot(3,1,1);
title(sprintf('Electrons',plotPop, popNum));
xlabel('x (nm)');
ylabel('y (nm)');
axis([0 L/1e-9 0 W/1e-9]);
hold on;
for i=1:plotPop
    plot(traj(:,i*2)./1e-9, traj(:,i*2+1)./1e-9, '.');
    
end


for j=1:size(boxes,1)
   plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)]./1e-9,...
       [boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)]./1e-9, 'k-');
end

subplot(3,1,3);
v = sqrt(state(:,3).^2 + state(:,4).^2);
title('Histogram of Electron Speeds');
histogram(v);
xlabel('Speed (m/s)');
ylabel('Number of particles');



density = hist3(state(:,1:2),[200 100])';

N = 20;
sigma = 3;
[x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(4);
imagesc(conv2(density,f,'same'));
set(gca,'YDir','normal');
title('Electron Density');
xlabel('x (nm)');
ylabel('y (nm)');


tempSumX = zeros(ceil(L/1e-9),ceil(W/1e-9));
tempSumY = zeros(ceil(L/1e-9),ceil(W/1e-9));
tempSum = zeros(ceil(L/1e-9),ceil(W/1e-9));


for i=1:popNum
   
    x = floor(state(i,1)/1e-9);
    y = floor(state(i,2)/1e-9);
    if(x==0)
        x = 1;
    end
    if(y==0)
        y= 1;
    end
    
  
    tempSumY(x,y) = tempSumY(x,y) + state(i,3)^2;
    tempSumX(x,y) = tempSumX(x,y) + state(i,4)^2;
    tempSum(x,y) = tempSum(x,y) + 1;
end



temp = (tempSumX + tempSumY).*m./k./2./tempSum;
temp(isnan(temp)) = 0;
temp = temp';



N = 20;
sigma = 3;
[x y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
f=f./sum(f(:));
figure(5);
imagesc(conv2(temp,f,'same'));
set(gca,'YDir','normal');
title('Temperature Map');
xlabel('x (nm)');
ylabel('y (nm)');




