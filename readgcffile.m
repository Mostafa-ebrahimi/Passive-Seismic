function [samples,streamID,sps,ist] = readgcffile(filename, streamID)
% ReadGCFFile
%
%   [SAMPLES,STREAMID,SPS,IST] = READGCFFILE(filename, streamID)
%
%   Reads in the specified GCF formatted file, and returns:
%     Samples - an array of all samples in file
%     Stream ID (string up to 6 characters)
%     SPS - sample rate of data in SAMPLES
%     IST - start time of data, as serial date number
%
%   example:
%   [samples,streamID,sps,ist]=readgcffile('test.gcf');
%   streams=readgcffile('test.gcf','list');
%   [samples,streamID,sps,ist]=readgcffile('test.gcf','TESTZ2');
%
%   M. McGowan, Guralp Systems Ltd.
%   2004/09/23 M. McGowan (support@guralp.com)
%     Added support for multiple streams in a GCF file, where the user can
%     specify which stream ID they want to extract. See 'streamID' input
%     parameter. This is optional - if ommited, it will use the first
%     streamID it finds. Note that it IS case-sensitive - all IDs should be
%     uppercase. If StreamID is a cell array of strings, this function will
%     return array structures containing data for all streams specified.
%     [SAMPLES,STREAMID,SPS,IST] = READGCFFILE(FILENAME, 'list')
%     Specifying 'list' for the streamID will return a cell array of strings,
%     one string for each streamID found in the file. This can be used to
%     iterate through each streamID in the file to read all the data contained
%     in the file. For an example, see the 'plot' option.
%     READGCFFILE(FILENAME, 'plot')
%     The 'plot' option is an example of reading all streams in a file and
%     displaying them.
%
%     If the specified stream is a status stream, 'samples' will return an
%     array of numbers which can be converted into text using
%     char(samples').
%
%     Modified code to cope with gaps, overlaps and out-of-sequence data.
%     Uses the first block timestamp as a reference, so will not return any
%     data in the file that has a timestamp older than the first
%     block in file.
%     In the case of an overlap, the data found further through the file will
%     overwrite the data read earlier.
%     Where a gap exists, it will be padded with sample values of NaN.
%
%   2008/11/13 M.McGowan (support@guralp.com)
%     Updated plotfile routine - only plots data streams (sps>0), and plots
%     aligned for time
%
%   2008/11/28 M.McGowan (support@guralp.com)
%     Updated readgcfblock for sample rates >250
%
%   2009/03/05 M.McGowan (support@guralp.com)
%     Fixed bug picking up streamID from first block. Was preventing code
%     from reading a GCF file without a Stream ID being specified


% initialise output variables to 'invalid'. Can't leave them non-existing,
% as that doesn't compile properly in 6.x.
sps=-1;
ist=-1;
wstreamID=-1;
expectedtime=-1;

if (nargin>1) && iscell(streamID),
  for i = 1:length(streamID),
    [samples{i},streamID{i},sps(i),ist(i)]=readgcffile(filename,streamID{i});
  end
  return
end

if (nargin>1) && strcmp(streamID,'plot'),
  plotfile(filename);
  return
end

fid = fopen(filename,'r','ieee-be');
if fid==-1,
    [p,n,e]=fileparts(filename);
    if ~strcmpi(e,'.gcf'),
        fname2 = [filename,'.gcf'];
        fid = fopen(fname2,'r','ieee-be');
    end
end
if fid==-1,
    error(['Unable to open file "',filename,'"']);
    return; 
end

if (nargin>1) && strcmp(streamID,'list'),
  samples=getstreamidlist(fid);
  fclose(fid);
  return
end

if nargin>1,
  wstreamID = base2dec(streamID,36); % faster to compare numbers than strings, so use a 'working' streamID
end
% to read the file, first create the array to handle the entire file's samples,
% then read in block by block, copying into the array in the correct place.
% This is MUCH faster than adding on to the end of an array each block.
sampcount=samplesinfile(fid);
%samples=NaN*ones(sampcount,1);
samples=repmat(NaN,sampcount,1); % faster than the above, according to docs for repmat
sampcount=1;

onesec = datenum(0,0,0,0,0,1);
onemsec = onesec/1000;

while ~feof(fid)
  if wstreamID < 0, % no stream ID has been pre-specified, so use the first one we find
    [blksamples,blksysID,blkstreamID,blksps,blkist] = readgcfblock(fid);
    wstreamID = blkstreamID;
  else
    [blksamples,blksysID,blkstreamID,blksps,blkist] = readgcfblock(fid,wstreamID);
  end
  if sps < 0,
    sps = blksps;
  end
  if ist < 0,
    ist = blkist;
  end
  if (sps>0) && (~isempty(blksamples)),
    if expectedtime>=0,
      if expectedtime+onemsec < blkist,
        disp(['Warning: Gap in ',dec2base(blkstreamID,36),', Expected ',datestr(expectedtime,31),', found ',datestr(blkist,31)]);
      end
      if blkist+onemsec < expectedtime,
        disp(['Warning: Overlap in ',dec2base(blkstreamID,36),', Expected ',datestr(expectedtime,31),', found ',datestr(blkist,31)]);
      end
    end
    secs = length(blksamples)/sps;
    expectedtime = blkist + secs*onesec;
  end
  % Copy the samples into the pre-prepared array
  if blksps>0,
    ofs = round((blkist-ist)*blksps/onesec);
  else
    ofs = sampcount;
  end
  endofs = ofs + length(blksamples);
  while endofs > length(samples), %if array not big enough, expand until it is
    samples = [samples;NaN*ones(length(samples),1)];
  end
  
  if ofs>=0,
    samples(ofs+1:endofs)=blksamples;
    sampcount = max([sampcount,endofs]);
  else
    disp(['Warning: discarding data from',dec2base(blkstreamID,36),' as it is before the start of file. FileStart=',datestr(ist,31),' BlockStart=',datestr(blkist,31)]);
  end
end
fclose(fid);
samples=samples(1:sampcount);     % trim samples array to actual length
streamID = dec2base(wstreamID,36); % convert numerical streamID back to a string



function [samps,sysID,streamID,sps,ist] = readgcfblock(fid,nstrid)
samps=[];
sps=0;
ist=0;
sysID = fread(fid,1,'uint32');
streamID = fread(fid,1,'uint32');
if nargin>1, % if we have specified a particular ID, keep searching until we find it
  while ~feof(fid) && (nstrid ~= streamID),
    fseek(fid,1016,'cof');
    sysID = fread(fid,1,'uint32');
    streamID = fread(fid,1,'uint32');
  end
end
if feof(fid)
  return 
end

date = fread(fid,1,'ubit15');
time = fread(fid,1,'ubit17');
reserved = fread(fid,1,'uint8');
sps = decodesps(fread(fid,1,'uint8'));
frac = fread(fid,1,'ubit4');
compressioncode = fread(fid,1,'ubit4');
numrecords = fread(fid,1,'uint8');

% Convert GCF coded time to Matlab coded time
hours = floor(time / 3600);
mins = rem(time,3600);
ist = datenum(1989,11,17, hours, floor(mins / 60), rem(mins,60) ) + date;
% add in the fractional second offset (if any)
if sps>0,
  if sps==400, step=50; else step=250; end;
  diff=frac*step/sps; % fractions of a second
  ist=ist+ diff/86400;
end

if (sps ~= 0),
   fic = fread(fid,1,'int32');
   switch compressioncode
   case 1,
      diffs = fread(fid,numrecords,'int32');
   case 2,
      diffs = fread(fid,numrecords*2,'int16');
   case 4,
      diffs = fread(fid,numrecords*4,'int8');
   end
   ric = fread(fid,1,'int32',1000-numrecords*4);
   diffs(1) = fic;
   samps = cumsum(diffs);
else
   samps = char(fread(fid,numrecords*4,[num2str(numrecords*4),'*uchar=>uchar'],1008-numrecords*4)');
end




function samps = samplesinfile(fid)
fseek(fid,14,'bof');
% Read number-of-records and compression-code of every block into an array
nr = fread(fid,'uint16',1022);
% Separate number-of-records and compression-code from the 16 bit value read
cc = bitshift(nr,-8);
nr = bitand(nr,255);
% sum up the number of samples in each block
samps=sum(cc.*nr);
frewind(fid);



function list = getstreamidlist(fid)
fseek(fid,4,'bof');
list=fread(fid,'uint32',1020);
list=dec2base(list,36);
list=unique(cellstr(list));


function plotfile(fname)
% EXAMPLE SCRIPT TO READ AND PLOT ALL STREAMS IN A GCF FILE
list=readgcffile(fname,'list');
for i = 1:length(list),
  [samples{i},id{i},sps(i),Tstart(i)]=readgcffile(fname,list{i});
  Tend(i)=Tstart(i);
  if sps(i)>0,
    Tend(i)=Tstart(i) + (length(samples{i})/sps(i))/86400;
  end
end
Tmin=datevec(min(Tstart));
Tmax=datevec(max(Tend));
Tdiff=etime(Tmax,Tmin);
datachans=find(sps>0);
numchans=length(datachans);
for i = 1:numchans,
  subplot(numchans,1,i);
  Soffset=etime(datevec(Tstart(datachans(i))),Tmin);
  plot((Soffset:1/sps(datachans(i)):Soffset+(length(samples{datachans(i)})-1)/sps(datachans(i))),samples{datachans(i)});
  ylabel(id(datachans(i)));
  ax=axis;
  axis([0 Tdiff ax(3) ax(4)]);
end


function [outsps]=decodesps(insps)
switch insps
  case 157,  outsps=0.1;
  case 161,  outsps=0.125;
  case 162,  outsps=0.2;
  case 164,  outsps=0.25;
  case 167,  outsps=0.5;
  case 171,  outsps=400;
  case 174,  outsps=500;
  case 176,  outsps=1000;
  case 179,  outsps=2000;
  case 181,  outsps=4000;
  otherwise outsps=insps;
end