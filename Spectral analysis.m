%This code is wrriten by mostfa ebrahimi
%the master student of Geophysics, in university of Tehran
%This code is about V/H in Passive Seismic Studies
%----------------------------------------------------------------------
clc
clear 
close all
%----------------------------------------------------------------------
%mul is about changing fluctuation means... 
%this changes smoothing of V/H
%curve
mul=30;
%please enter the path which save the data
 cd('E:\sssssss\');
 %path the 3-component data
Zdata=dir('E:\sssssss\*z.gcf');
Edata=dir('E:\sssssss\*e.gcf');
Ndata=dir('E:\sssssss\*n.gcf');
%This is about gcffile for guralp seismometer
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
  %___________________________________________________________________________________________________    
    
  for ni=1:3
      if ni==1
          timein=sz;
      end
          if ni==2
              timein=sn;
          end
              if ni==3
                  timein=se;
              end
                  %--------------------------------------------------------------------------
                  cont=-1;
                  sps=100;
                  CMG3ESP.p = [(-23.65e-3)+23.65e-3i;(-23.65e-3)-2365e-3i;-180;-160;-80] * 2 * pi;          %rad/s
                  CMG3ESP.z = [0 0];                                                                     %rad/s

                  nf_G3ESP = 2304000;            %At 1 Hz

                  ss_G3ESP = 2 * 2979;            %volt.sec/m (velocity output)

                  sd_G3ESP = 3.20e-6;      %volt/count

                  CMG3ESP.k = nf_G3ESP * ss_G3ESP / sd_G3ESP * (2 * pi)^3             %count.sec/m (for CMG-DM24S3)
                  %--------------------------------------------------------------------------

                  %-----------------------------------------------------------------------

                  

                  len = 2.^nextpow2(length(timein));

                  %Create a frequency vector
                  fvec = sps * linspace(0,1,len)';

                  %Change to the angle frequency
                  omeg = zeros(size(fvec)) + 1i * 2 * pi * fvec;

                  %Pick out the poles and the zeros and the amplitude
                  %amp = resp.const;
                  amp=CMG3ESP.k;
                  %zerores = resp.zer;
                  zerores=CMG3ESP.z;
                  %poleres = resp.pol;
                  poleres=CMG3ESP.p;

                  %Calculate the zeros-response for different frequency
                  for k = 1:length(zerores)
                      tem0 = real(zerores(k)) * ones(len,1) + 1i * imag(zerores(k)) * ones(len,1);
                      tem1(:,k) = omeg - tem0;
                  end
                  tem1 = tem1';
                  zinst = prod(tem1)';

                  clear tem1 tem0

                  %Calculate the poles-response for different frequency
                  for k = 1:length(poleres)
                      tem0 = real(poleres(k)) * ones(len,1) + 1i * imag(poleres(k)) * ones(len,1);
                      tem1(:,k) = omeg - tem0;
                  end
                  tem1 = tem1';
                  pinst = prod(tem1)';

                  clear tem1 tem0

                  %zeros and poles response are put together

                  norm2in = abs(pinst) .^ 2;
                  maxv = max(max(norm2in));

                  %Determine the waterlevel according to the specific case
                  temp = norm2in;
                  whozero = find(temp == 0);
                  temp(whozero) = maxv * ones(size(whozero));
                  minv = min(min(temp));
                  norm2in(whozero) = minv * ones(size(whozero));

                  instres = amp * (zinst .* conj(pinst)) ./ norm2in;

                  clear norm2in

                  %figure,loglog(fvec,abs(instres));

                  %Work with data
                  if cont == 1 % Add the instrument response
                      fftdat = fft(timein,len);
                      fftime = fftdat .* instres;
                  elseif cont == -1 % Remove the instrument response
                      fftdat = fft(timein,len);
                      norm2in = abs(instres).^2;
                      maxv = max(max(norm2in));

                      %Determine the waterlevel according to the specific case
                      temp = norm2in;
                      whozero = find(temp == 0);
                      temp(whozero) = maxv * ones(size(whozero));
                      minv = min(min(temp));
                      norm2in(whozero) = minv .* ones(size(whozero));

                      fftime = (fftdat .* conj(instres)) ./ norm2in;
                  end

                  clear norm2in fftdat

                  % Return to the time domain
                  timeout = real(ifft(fftime));

                  %Original length of the time series is recovered
                  timeout = timeout(1:length(timein));
                  if ni==3
                  se2=timeout;
                  end
              if ni==2
              sn2=timeout;
              end
          if ni==1
          sz2=timeout;
          end

   end
  

%__________________________________________________________________________

%     snn=((sn).^2+(se).^2).^.5;
      %---------------------
    Fs=sps; Wn=[1 50]/(Fs);n=2;
    [b,a]=butter(n,Wn);
    velz=filter(b,a,sz2);
    vele=filter(b,a,se2);
    veln=filter(b,a,sn2);
%     velnn=filter(b,a,snn);
    [fz, magz, phasez] = spect(velz, Fs);
    [fe, mage, phasee] = spect(vele, Fs);
    [fn, magn, phasen] = spect(veln, Fs);
%     [fnn, magnn, phasenn] = spect(velnn, Fs);
    rec{i}.fz=fz;
    rec{i}.magz=magz;
    rec{i}.mage=mage;
    rec{i}.magn=magn;
%     rec{i}.magnn=magnn;
    
end
% spectview(velz)
cd ..
VH=[];
sm=mul*Fs;
%--------------------------------------------------------------------------
for p=1:length(rec)
    spZ=smooth(rec{i}.magz,sm);
    spE=smooth(rec{i}.mage,sm);
    spN=smooth(rec{i}.magn,sm);
%     spnN=smooth(rec{i}.magnn,sm);
      voh=1./(0.5*((spN+spE))./spZ);
%  voh=spZ./(spE);
%     voh=1./((spnN)./spZ);
    voh=smooth(voh,sm);
    fL= rec{i}.fz;
%     figure(1)
%     plot(fL,voh);
%     hold on
%     
%       figure(3)
%     plot(fL,spE);
%     hold on
    
%     xlim([0.01,Fs/2]);
    VH=[VH,voh];
end
vhmean=mean(VH,2);
% index = 1:4;
% baseLine=0.2;
% % area(fL<4 & fL>1 ,voh);
figure(2)
plot(fL,vhmean,'k', 'LineWidth',3);
%  fill(fL,vhmean,'r')
% hold on;  
% plot(fL,vhmean,'b', 'LineWidth',2);
%  patch(fL,vhmean,[1 1 0]);
% set(fL,'BaseValue',-2)
  hold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim([0.2,10]);
ylim([2.622,2.6245])
% xrange=2:4;
bb=fix(4/(max(fL)/length(fL)));
aa=fix(2/(max(fL)/length(fL)));
area(fL(aa:bb),vhmean(aa:bb));
% ymin=2.622;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot(x,y1,'b');                              %# Plot the first line
% hold on;                                     %# Add to the plot
% h1 = fill(fL(aa:bb),...        %# Plot the first filled polygon
%           [vhmean(aa:bb)],...
%           'g','EdgeColor','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
xlabel('Frequency (Hz)');
ylabel('V/H')

%  h1 = fill(fL(index([1 1:end end])),...        %# Plot the first filled polygon
%           [baseLine vhmean(index).' baseLine],...
%           'r','EdgeColor','none');



text(3.1765   , 2.62351,'*','fontsize',39);



% Create line
% for vhmeani=1:length(vhmean)
%     if ll==vhmean
%     ll=vhmean;
%     end
% end
% 
% [fL,vhmean]=size(a);
% for i=1:length(vhmean)
%     max=a(i,1)
%     for j=i:n
%         if a(i,j)>max
%             max=a(i,j)
%         end
%     end
%     m(i)=max
% end
    
    
    
% annotation(figure(2),'line',[0.4305 ll],[0.1101 0.4654],'LineWidth',2,...
%     'Color',[1 0 0]);
% 
% % Create line
% annotation(figure(2),'line',[0.194 0.194],[0.1101 0.381],'LineWidth',2,...
%     'Color',[1 0 0]);
%Area
% h2 = area(fL,vhmean,'FaceColor',[.5 .9 .6],'EdgeColor','b', 'LineWidth',2);
% xx = fill(x(index([1 1:end end])), [baseLine y1(index) baseLine],'b','EdgeColor','none');

% y=vhmean
% figure(3)
% 
% fL=Linspace(1,10,262144);
% area(fL,y)
% patch(fL,y,[1 1 0])

%   ahmad

% createfigure(fL1, VH1);

sddd = sum (vhmean (26000:52000));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hvmean = 1./vhmean;
% figure, plot (fL,hvmean)
% figure, plot(fL,hvmean,'b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim([0.2,10]);
ggg = vhmean (44564:96994);
ggg2 = ggg - 1;
ggg3 = ggg2 (ggg2 > 0);

% xs=fL(fL>1 & fL<4);
% figure;
% hold on;
% area(xs,3);
% % area(xs,normpdf(xs,0,3));
% plot(x,3);


M = vhmean(26214);

% for i=1:length(vhmean)-1
for i=26214:36993
    
 if vhmean(i+1)>=M
        M = vhmean(i+1);
 end
    
end
disp('max v/h=')
disp(M)

% bar(fl,vhmean)
