%This code is wrriten by mostfa ebrahimi
%the master student of Geophysics, in university of Tehran
%This code is about V/H in Passive Seismic Studies
%----------------------------------------------------------------------
clc
clear all
sps=100;
cd('D:\');
Zdata=dir('D:\*Z_MSEED');
Ndata=dir('D:\*N_MSEED');
Edata=dir('D:\*E_MSEED');

 for k=1:length(Zdata)
    clear  spZ spN spE l1 l 
    kz=Zdata(k).name;
    kn=Ndata(k).name;
    ke=Edata(k).name;

z=rdmseed(kz);
n=rdmseed(kn);
e=rdmseed(ke);
[z1 z2]=size(z);
[n1 n2]=size(n);
[e1 e2]=size(e);
sz=[];
sn=[];
se=[];
for i=1:z2
    sz=[sz;z(1,i).d];
end

for j=1:n2
    sn=[sn;n(1,j).d];
 end

for q=1:e2
    se=[se;e(1,q).d];
end

 l1=10000;  
         l=size(sz)/l1-1;
    
    dcn=mean(sn);
    dce=mean(se);
    dcz=mean(sz);
    sn=sn-dcn;
    se=se-dce;
    sz=sz-dcz;
     
     for j3=1:floor(l(:,1));
%          %---------------------
        Fs=sps; Wn=[0.5 12]/(100);n=4;
        [b,a]=butter(n,Wn);
        
        dcn=mean(sn((l1*(j3-1)+1):(j3*l1)));
        dce=mean(se((l1*(j3-1)+1):(j3*l1)));
        dcz=mean(sz((l1*(j3-1)+1):(j3*l1)));
        sn2=sn-dcn;
        se2=se-dce;
        sz2=sz-dcz;

     end
     

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

                  nf_G3ESP = 2304000;                                                               %At 1 Hz

                  ss_G3ESP = 2 * 2979;                                                              %volt.sec/m (velocity output)

                  sd_G3ESP = 3.20e-6;                                                               %volt/count

                  CMG3ESP.k = nf_G3ESP * ss_G3ESP / sd_G3ESP * (2 * pi)^3                          %count.sec/m (for CMG-DM24S3)

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
  

 end
 
 
 g=sn;
plot(g,'b')

tseries=sn;dt=0.01;stw=2;ltw=100*stw;thresh=8;
ln = fix( ltw / dt);
shn = fix(stw / dt);
nt = length(tseries);
aseries = abs( hilbert(tseries) );
sra = zeros(1, nt);
for aa = ln + 1:nt
    lta = mean(aseries(aa - ln : aa));
    sta = mean(aseries(aa - shn : aa));
    sra1(aa) = sta / lta;
end
figure;plot(sn) 
figure;plot(sra1) 
itm = find(sra > thresh);
if ~isempty(itm)
    itmax = itm(1);
end
%-----------------------------------------------------------------------
tseries=sz;dt=dcz;stw=100*dt;ltw=10*stw;thresh=5;
ln = fix( ltw / dt);
shn = fix(stw / dt);
nt = length(tseries);
aseries = abs( hilbert(tseries) );
sra = zeros(1, nt);
for aa = ln + 1:nt
    lta = mean(aseries(aa - ln : aa));
    sta = mean(aseries(aa - shn : aa));
    sra2(aa) = sta / lta;
end
figure;plot(sz) 
figure;plot(sra2) 
itm = find(sra > thresh);
if ~isempty(itm)
    itmax = itm(1);
end
%-----------------------------------------------------------------------
tseries=se;dt=dce;stw=100*dt;ltw=10*stw;thresh=5;
ln = fix( ltw / dt);
shn = fix(stw / dt);
nt = length(tseries);
aseries = abs( hilbert(tseries) );
sra = zeros(1, nt);
for aa = ln + 1:nt
    lta = mean(aseries(aa - ln : aa));
    sta = mean(aseries(aa - shn : aa));
    sra3(aa) = sta / lta;
end
figure;plot(se) 
figure;plot(sra3) 
itm = find(sra > thresh);
if ~isempty(itm)
    itmax = itm(1);
end

S=[se';sn';sz'];
SRA=[sra1;sra2;sra3];
figure;imagesc(S)
figure;imagesc(SRA)

