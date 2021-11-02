%Matlab code for simulating the trajectories of bacteria.  

clear all
clf

%parameters
Dr=0.5; %rotational diffusivity
flip=3; %reversal frequency
N=100; %number of particles
U=50; %speed in microns/second
T=10; %total time in seconds
F=T*25; %number of frames
delay=12;  %minimum number of frames between successive reversals

dt=T/F;  %time step
time=0:dt:T-dt;    

%generation of random numbers
b1=randn(N,F);
b2a=random('Geometric',dt*flip,N,1)+1;      %first reversal without delay
b2b=random('Geometric',dt*flip,N,F)+delay;    %subsequent reversals have delay
b2=zeros(N,F);

for i=1:N
    
    ii=b2a(i);
    b2(i,ii)=1;

    jj=1;
    while ii+b2b(i,jj)<F
        
        ii=ii+b2b(i,jj);
        b2(i,ii)=1;
        jj=jj+1;
    end
end 

%Simulated data
posx=zeros(N,F);   %x coordinates
posy=zeros(N,F);   %y coordinates 
ang=zeros(N,F);    % heading angles 

%initiate all bacteria to enter the channel at x=0
posx(:,1)=zeros(N,1);
posy(:,1)=zeros(N,1);
ang(:,1)=zeros(N,1);

ind=1:N;  %indices of particles remaining in the channel x>0

j=1;
while j<F & length(ind)>0
    
    %time-marching step
    j=j+1;
             
    posx(ind,j)=posx(ind,j-1)+U*cos(ang(ind,j-1))*dt;
    posy(ind,j)=posy(ind,j-1)+U*sin(ang(ind,j-1))*dt;
    ang(ind,j)=ang(ind,j-1)+pi*b2(ind,j-1)+sqrt(2*Dr*dt)*b1(ind,j-1);
    
    p1=find(posx(ind,j)<0);
    %bacteria leaving the channel are identified
    if length(p1)>0
        
        residence(ind(p1))=j;
        posx(ind(p1),j+1:F)=repmat(posx(ind(p1),j),1,F-j);
        posy(ind(p1),j+1:F)=repmat(posy(ind(p1),j),1,F-j);
        ind(p1)=[];
    end
    
    p2=find(abs(posy(ind,j))>50);
    %bacteria crossing the side walls are set to remain in the channel
    posy(ind(p2),j)=50*sign(posy(ind(p2),j));
    ang(ind(p2),j)=sign(cos(ang(ind(p2),j)))*pi/2-pi/2;   
            
end

%Plotting trajectories and their mean

plot(time,mean(posx,1),'r-','Linewidth',2)
hold on

for i=1:N
    
    plot(time(2:end),posx(i,2:end),'k-')
    
end

xlabel('Time (s)')
ylabel('x-distance (um)')
set(gca,'Fontsize',20)
axis([0 10 0 150])
box on


