%This code is wrriten by Mostfa Ebrahimi
% Master student of Geophysics,  university of Tehran
%This code is about Polarization for passive Seismic
%----------------------------------------------------------------------
clc
close all
clear all
sps=100;
VH=[];
vohs=[];
% cd('D:\thesis\mar5\1108-1116\');
Zdata=dir('*z.gcf');
Edata=dir('*e.gcf');
Ndata=dir('*n.gcf');
for i=1:length(Zdata)
    [zz,ID,sps,ist] = readgcffile(Zdata(i).name);
    sz=zz;
%     [zz(140000:268000)];
    [nn,ID,sps,ist] = readgcffile(Ndata(i).name);
    [ee,ID,sps,ist] = readgcffile(Edata(i).name);
    sn=nn;
%     [nn(140000:268000)];
    se=ee;
    
    
    dcn=mean(sn);
    dce=mean(se);
    dcz=mean(sz);
    sn=sn-dcn;
    se=se-dce;
    sz=sz-dcz;
%     [ee(140000:268000)];
  %_
        sn2=sn-dcn;
        se2=se-dce;
        sz2=sz-dcz;
     end

        sz2=sz2(70000:100000);
         se2=se2(70000:100000);
          sn2=sn2(70000:100000);
% % % % % % % % % % pol
% % % % % % % % % % pol
t=0:0.01:(size(sz2)/100)-0.01;
t=t';
z=sz2;
n=sn2;
e=se2;
flp=0.01;
fhi=0.07;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%figure 1%%%%%%
if flp > 0;
   w = [flp fhi];
   [b,a]=butter(4,w);
 f2=figure('name','Filterd seismograms'); 
 % plot the filtered data 
 subplot(3,1,1); 
 plot(t,e); 
 xlabel('time (s)');
%  ylabel(strcat('EW Comp. of ',IDe)); 
 subplot(3,1,2); 
 plot(t,n); 
 xlabel('time (s)'); 
%  ylabel(strcat('NS Comp. of ',IDn)); 
 subplot(3,1,3); 
 plot(t,z); 
 xlabel('time (s)'); 
%  ylabel(strcat('Z comp. of ',IDz)); 

else 
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Moving window loop 
delt=1/sps;
ttot=600;
twin=1;
npts=ttot/delt + 1; 
npts1=fix(((size(z)-1)*1/sps)/delt) + 1;  
%  total number of samples to analyze 
nwin=fix(twin/delt) + 1 ;  
%  number of samples in a time window 
npshift=fix(twin/(2*delt))+1 ;
% number of samples to shift over 
kfin=fix((npts1-nwin)/(npshift+1))+1;
% number of time windows considered 
mxde1=0; 
mxde2=0; 
mxde3=0; 
for k=1:kfin; 
   nwinst=(k-1)*(npshift-1)+1;  
   % start of time window 
   nwinfn=nwinst+nwin-1;   
   % end of time window 
   a=[];
   a(:,1)=e(nwinst:nwinfn);
   a(:,2)=n(nwinst:nwinfn);
   a(:,3)=z(nwinst:nwinfn);   
   %a 
   %  c = cov(a(:,1),a(:,2),0);
   c=(1/length(e))*(a'*a);      
   % correlation matrix 
   %c 
   % c=a'*a;   % correlation matrix
   [v,d]=eig(c);      
   % eigen vectors and eigen values
   d ;    %eigen values
   v;
   %--------------------------------------------------%
    % sort the eigenvalues and eigenvectors 
   
    % Azimuth for first  eigenvalue
   ang1(k)=atan2(v(1,3),v(2,3)) * 180/pi;
   %_-fo---d-f--oflim652l>546
   ang2(k)=atan2(v(1,2),v(2,2)) * 180/pi;
   % Azimuth for second eigenvalue 
   ang3(k)=atan2(v(1,1),v(2,1)) * 180/pi; 
   % Azimuth for third  eigenvalue
   %---------------------------------------------------%
   dip1(k)=atan2(v(3,3),((v(1,3)^2+v(2,3)^2)^0.5)) * 180/pi;
   % Dip for first  component
   %---------------------------%
   dip2(k)=atan2(v(3,2),((v(1,2)^2+v(2,2)^2)^0.5)) * 180/pi;
   % Dip for second component 
   %----------------------------%
   dip3(k)=atan2(v(3,1),((v(1,1)^2+v(2,1)^2)^0.5)) * 180/pi;
   % Dip for third  component
   %------------------------------------------------------%
   de1(k)=d(1) ;
   de2(k)=d(2); 
   de3(k)=d(3); 
   LL(k)=1-((d(2,2)+d(1,1))/(d(3,3)*2));
   DD(k)=d(3,3);
   % find the maximum values 
   mxde1=max(mxde1,de1(k));   
   mxde2=max(mxde2,de2(k)); 
   mxde3=max(mxde3,de3(k)); 
   %angle from the vertical
 vang1(k)=acos(abs(v(3,1)))* 180/pi;   
 vang2(k)=acos(abs(v(3,2)))* 180/pi;  
 vang3(k)=acos(abs(v(3,3)))* 180/pi;
 t2(k)=delt*(nwinst-1);  
 % assign time for this window to the window start 
end;
L = 1 - (d(2,2)+d(3,3))/d(3,1)/2;
%-------------------------------------------------

%%%%%figure 2%%%%


f4=figure('name','Eigenvalues and Inferred Azimuth'); 
subplot(4,1,1); 
% plot(t2,de1,'-or',t2,de2,'-dg',t2,de3,'-+b');
plot(t2,DD);
xlabel('time sec'); 
ylabel('Strength');
%--------------------------------------------

subplot(4,1,2); 
%  plot(t2,ang1,'-or',t2,ang2,'-dg',t2,ang3,'-+b'); 
% plot(t2,ang3); 
plot(t2,ang1);
xlabel('time sec'); 
ylabel('Azimuth '); 
%-------------------------------------------------------%
subplot(4,1,3); 
%  plot(t2,vang1,'-or',t2,vang2,'-dg',t2,vang3,'-+b'); 
plot(t2,dip1);
xlabel('time sec'); 
ylabel('Dip ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,1,4); 
%  plot(t2,vang1,'-or',t2,vang2,'-dg',t2,vang3,'-+b'); 
plot(t2,LL);
xlabel('time sec'); 
ylabel('Rectilinearity');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalues and Inferred Azimuth
f3=figure('name','Eigenvalues and Inferred Azimuth'); 
subplot(2,1,1); 
% plot(t2,de1,'-or',t2,de2,'-dg',t2,de3,'-+b');
plot(t2,de1,'-+b');
xlabel('time sec'); 
ylabel('eigenvalues');

%-------------------------------------------------------%
subplot(2,1,2); 
%  plot(t2,vang1,'-or',t2,vang2,'-dg',t2,vang3,'-+b'); 
plot(t2,vang3,'-+b');
xlabel('time sec'); 
ylabel('incidence angle ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compass plots 
 f5=figure('name','Compass Plots'); 
%-------------------------------------------------%
subplot(2,3,4); 
compass(de1.*cos(ang1*pi/180),de1.*sin(ang1*pi/180)); 
title('Azimuth - Largest Eigenvalue','FontSize',8); 
%-----------------------------------------------%
subplot(2,3,5); 
compass(de2.*cos(ang2*pi/180),de2.*sin(ang2*pi/180)); 
title('Azimuth - Intermediate Eigenvalue','FontSize',8); 
%-----------------------------------------------%
subplot(2,3,6);  
compass(de3.*cos(ang3*pi/180),de3.*sin(ang3*pi/180));
title('Azimuth - Smallest Eigenvalue','FontSize',8);
%------------------------------------------------%
% eigenvalue 
nskip=1; 
if nskip == 1; 
neig1=0; 
neig2=0; 
neig3=0; 
for k=1:kfin; 
   if de1(k) >= 0.3*mxde1 && de1(k) < 0.999*mxde1; 
      neig1=neig1+1; 
      angm1(neig1)=ang1(k); 
   else 
   end; 
   if de2(k) >= 0.3*mxde2; 
      neig2=neig2+1; 
      angm2(neig2)=ang2(k); 
   else 
   end; 
 if de3(k) >= 0.3*mxde3; 
      neig3=neig3+1; 
      angm3(neig3)=ang3(k); 
   else 
   end; 
end; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rose plot 
f4=figure('name','Azimuth Distribution'); 
%-----------------------------------------------------%
subplot(2,3,1); 
rose(ang1*pi/180,100); 
title('Az- Largest Eig');
%-------------------------------------------------------%
subplot(2,3,2); 
rose(ang2*pi/180,100); 
title('Az- Intermediate Eig'); 
%------------------------------------------------------%
subplot(2,3,3); 
rose(ang3*pi/180,100); 
title('Az- Smallest Eig'); 
%------------------------------------------------------%
subplot(2,3,4); 
rose(angm1*pi/180,100); 
 title('Az-L.Eig,50% Threshold','FontSize',10);
%  title('Azimuth - Largest Eigenvalue,50% Threshold','FontSize',8);
%-----------------------------------------------------%
subplot(2,3,5); 
rose(angm2*pi/180,100);
title('Az-In.Eig,50% Threshold','FontSize',10);
% title('Azimuth - Intermediate Eigenvalue,50% Threshold','FontSize',8); 
%------------------------------------------------------%
subplot(2,3,6);  
rose(angm3*pi/180,100); 
title('Az-S.Eig,50% Threshold','FontSize',10);
% title('Azimuth - Smallest Eigenvalue,50% Threshold','FontSize',8);
%--------------------------------------------------------%