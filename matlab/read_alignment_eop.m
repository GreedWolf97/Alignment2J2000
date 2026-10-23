function q=read_alignment_eop(file,utc)
% Read CelesTrak EOP v1.1 observed and predicted data (13 columns).
fid=fopen(file,'r');
if fid<0, error('Alignment:EOPFile','Cannot open %s',file); end
clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
rows=zeros(0,13); pred=false(0,1); section='';
while true
    line=fgetl(fid); if ~ischar(line), break; end
    line=strtrim(line);
    if strcmp(line,'BEGIN OBSERVED'), section='observed'; continue; end
    if strcmp(line,'BEGIN PREDICTED'), section='predicted'; continue; end
    if startsWith(line,'END '), section=''; continue; end
    if isempty(section), continue; end
    v=sscanf(line,'%f');
    if numel(v)~=13 || any(~isfinite(v))
        error('Alignment:EOPFormat','Invalid EOP row: %s',line);
    end
    rows(end+1,:)=v.'; pred(end+1,1)=strcmp(section,'predicted'); %#ok<AGROW>
end
if isempty(rows)||any(diff(rows(:,4))<=0)
    error('Alignment:EOPFormat','EOP MJD must be strictly increasing.');
end
m0=juliandate(utc(1),utc(2),utc(3))-2400000.5;
f=(utc(4)*3600+utc(5)*60+utc(6))/86400; m=m0+f;
if m<rows(1,4)||m>rows(end,4)
    error('Alignment:EOPRange','MJD %.9f outside [%.0f,%.0f]; no extrapolation.',m,rows(1,4),rows(end,4));
end
k=find(rows(:,4)<=m,1,'last'); j=k;
if m>rows(k,4), j=k+1; end
if rows(k,4)~=m0 || (j~=k && rows(j,4)~=m0+1)
    error('Alignment:EOPGap','Missing daily EOP record.');
end
w=0; if j~=k, w=f; end
v=(1-w)*rows(k,5:12)+w*rows(j,5:12); dat=rows(k,13);
% DAT is a step. Interpolate UT1-TAI across a leap-second boundary.
du=(1-w)*(rows(k,7)-dat)+w*(rows(j,7)-rows(j,13))+dat;
q=struct('MJD',m,'xp',v(1),'yp',v(2),'deltaUT1',du,'LOD',v(4), ...
    'dPsi',v(5),'dEpsilon',v(6),'dX',v(7),'dY',v(8),'deltaAT',dat, ...
    'IsPredicted',pred(k)||(w>0&&pred(j)),'File',char(file));
end
