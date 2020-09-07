%%%%%%%%%% HILL'S VORTEX STREAKS GENERATION CODE %%%%%%%%%%
function I_pm =HillsVortexStreaksGeneration (PASettings)
%   PARAMETER'S SETTINGS
PASettings.Particle_Number = 10000;
PASettings.meanI = 0.5;  % (Normalized Intensity between 0~1.0)
PASettings.stdI = 0.1;
PASettings.meanR = 1; % (Radius)
PASettings.stdR = 0.25;
PASettings.I_size = [1000 1000];
PASettings.x_pixels_per_mm = 1; % (Magnification-x)
PASettings.y_pixels_per_mm = 1; % (Magnification-y)
ParticleList=Vortex (PASettings);
I_size=PASettings.I_size;
x_pixels_per_mm=PASettings.x_pixels_per_mm;
y_pixels_per_mm=PASettings.y_pixels_per_mm;
%   INITIALIZE
I_pm=zeros (I_size (1), I_size (2));
%   GENERATE PARTICLE IMAGE FROM PARTICLE LIST
for i=1:size(ParticleList,2)
   	Amp=ParticleList(i). Amp;
    	stdM=x_pixels_per_mm*(ParticleList(i). stdM);
    	stdm=y_pixels_per_mm*(ParticleList(i). stdm);
%   DETERMINE THE SIZE OF THE PARTICLE IMAGE
   	r_width=2*ceil(48*stdM) +1;
    	c_width=2*ceil(48*stdM) +1;
   	I_pm_local_size=[size(I_pm,1)+r_width-1 	size(I_pm,2)+c_width-1];
    	I_pm_local=zeros(I_pm_local_size);    
%   SHIFT THE COORDINATES BASED ON IMAGE COORDINATES (ORIGIN AT TOP-LEFT CORNER  
%   OF ZERO-PADDED I_PM LOCAL)
		 xc=ParticleList(i).px+(r_width-1)/2;
    	 yc=ParticleList(i).py+(c_width-1)/2;    
%   CALCULATE WHICH PIXEL THE PEAK CENTER FALLS INTO
   	 rc=round(xc)+1;
    	 cc=round(yc)+1;    
%   CALCULATE CENTER COORDINATES FOR MASK IMAGE    
    	 xc_mask=xc-(rc);
   	 yc_mask=yc-(cc);
    	py=ParticleList(i).px;
    	px=ParticleList(i).py;
%  VELOCITY USED TO GENERATE HILL'S VORTEX
     U=15;  % (pixel/interval)      
%  RADIUS USED TO GENERATE HILL'S VORTEX    
     R=700; % (pixel) 
%  Vx-Vy VELOCITY COMPONENTS OF THE HILL'S VORTEX
     Vx= U*((px-500).*(1000./500./R)).*((py-500).*(1000./500)./R);                 
     Vy= U*(1-((py-500).*(1000./500)./R).^2-2*((px-500).*(1000./500)./R).^2);       
     I_pm_local(rc-(r_width-1)/2:rc+(r_width-1)/2, cc-(c_width-1)/2:cc+ (c_width-1)/2) =ParticleMaskE ([r_width, c_width], Amp, stdM, stdm, Vx, Vy, xc_mask, yc_mask);
     I_pm=I_pm+I_pm_local((r_width-1)/2+1: I_pm_local_size (1) -(r_width-1)/2, (c_width-1)/2+1: I_pm_local_size (2) -(c_width-1)/2);
end
%   DISPLAY IMAGE WITH GAUSSIAN NOISE 5% OF THE MAXIMUM INTENSITY 
      figure (1); imagesc(I_pm);
   I_pm=imnoise(I_pm,'gaussian',0,0.0025);
%   PARTICLE MASK GENERATION FUNCTION
function [I_pm]=ParticleMaskE(I_size,Amp,stdM,stdm,Vx,Vy,varargin)
I_pm=ones(I_size(1),I_size(2));
[X,Y]=meshgrid(1:I_size(1),1:I_size(2));
if nargin>7
    xc=varargin{1}+double((I_size(1)+1)/2);
    yc=varargin{2}+double((I_size(2)+1)/2);
else
    xc=double((I_size(1)+1)/2);
    yc=double((I_size(2)+1)/2);
end
stdM=double(stdM);
stdm=double(stdm);
X=double(X);
Y=double(Y); 
%   u-v COMPONENT OF PARTICLE 
u=(X-xc)/stdM;
v= (Y-yc)/stdm;
% STREAK GENERATION BY TIME-INTEGRAL 2D GAUSSIAN FUNCTION
fun= @(t) Amp.*exp(-((u-Vx.*t).^2+(v-Vy.*t).^2)/2);
%   TIME-INTEGRAL (t)
t=2; 
G= integral(fun,0,t,'ArrayValued',true);
G=double(G);
I_pm=I_pm.*G;
%   FUNCTION TO GENERATE PARTICLE LIST AND HILL'S VORTEX 
function ParticleList=Vortex (PASettings)
N=PASettings.Particle_Number;
meanI=PASettings.meanI;
stdI=PASettings.stdI;
meanR=PASettings.meanR;
stdR=PASettings.stdR;
I_size=PASettings.I_size;
x_per_mm=PASettings.x_pixels_per_mm;
y_per_mm=PASettings.y_pixels_per_mm;
Pos=rand(N,2);
Radius=meanR+stdR*(rand(N,1)-0.5);
Amplitude=meanI+stdI*(rand(N,1)-0.5);
for i=1:N
    ParticleList(i).stdm=Radius(i);
    ParticleList(i).stdM=Radius(i);
    ParticleList(i).Amp=Amplitude(i);
    ParticleList(i).px=(Pos(i,1)*I_size(1))*x_per_mm;
    ParticleList(i).py=(Pos(i,2)*I_size(2))*y_per_mm;
    py=ParticleList(i).px;
    px=ParticleList(i).py;
    %   VELOCITY USED TO GENERATE HILL'S VORTEX   
    U=15;  % (pixel/interval)      
    %   RADIUS USED TO GENERATE HILL'S VORTEX
    R=700; % (pixel) 
    %   Vx-Vy VELOCITY COMPONENTS OF THE HILL'S VORTEX 
    Vx= U*((px-500).*(1000./500./R)).*((py-500).*(1000./500)./R);                
    Vy= U*(1-((py-500).*(1000./500)./R).^2-2*((px-500).*(1000./500)./R).^2);       
    ParticleList(i).py=px;
    ParticleList(i).px=py;
    ParticleList(i).Vx=Vx;
    ParticleList(i).Vy=Vy;
end 
