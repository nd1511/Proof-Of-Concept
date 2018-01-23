function [x]=interpofiltbankm_2(p,n,fs)
%INTERPOFILTBANKM determines interpolation matrix for a MEL to STFT bins transformation
%
%  VERY heavily inspired by FILTBANKM from Mike Brookes' voicebox MATLAB toolbox
%  (basically just changed the beginning to do the reverse operation, 
%   the rest of the code is pretty much the same)
%
% Inputs:
%       p   number of filters in filterbank or the filter spacing in k-mel/bark/erb [ceil(4.6*log10(fs))]
%		n   length of fft
%		fs  sample rate in Hz
%
% Outputs:	x     a matrix containing the filterbank amplitudes
%                 size(x)=[p,1+floor(n/2)]
%

%w='m';
w='b';
%wr ='m';
fh=0.5*fs; % max freq is the nyquist
fl=0;

f1=0;
nf=1+floor(n/2); % number of input frequency bins
df=fs/n;  % input frequency bin spacing
cf=f1+(0:nf)*df;  % input frequency bins

mflh=[fl fh];
%mflh=frq2mel(mflh);       % convert frequency limits into mel     
mflh=frq2bark(mflh);       % convert frequency limits into mel  
melrng=mflh*(-1:2:1)';          % mel/erb/... range
% fn2=floor(n/2);     % bin index of highest positive frequency (Nyquist if n is even)
melinc=melrng/(p+1);

fin0 = mflh(1)+(0:p+1)*melinc; % centre frequencies in mel/erb/... including dummy ends
fin0(2:end)=max(fin0(2:end),0); % only the first point can be negative
%fin0 = mel2frq(fin0);
fin0 = bark2frq(fin0);

cf = [cf(1)-df,cf];
%cf=mflh(1)+(0:p+1)*melinc; % centre frequencies in mel/erb/... including dummy ends
%cf(2:end)=max(cf(2:end),0); % only the first point can be negative

p = length(cf) - 2;

mb = cf;
%mb=bark2frq(cf);

% first sort out 2-sided input frequencies
fin=fin0;
%fin(end+1)=fin(end)+df; % add on a dummy point at the high end
fin=[-fin(end:-1:2) fin];
nfin=length(fin);  % length of extended input frequency list
nf = (nfin - 3)/2;

% now sort out the interleaving

fout=mb;  % output frequencies in Hz
lowex=any(w=='y')~=any(w=='Y');   % extend to 0 Hz
highex=any(w=='y') && (fout(end-1)<fin(end));  % extend at high end
if lowex
    fout=[0 0 fout(2:end)];
end
if highex
    fout=[fout(1:end-1) fin(end) fin(end)];
end
mfout=length(fout);
if any(w=='u') || any(w=='U')
    gout=fout(3:mfout)-fout(1:mfout-2);
    gout=2*(gout+(gout==0)).^(-1); % Gain of output triangles
else
    gout=ones(1,mfout-2);
end
if any(w=='u')
    gin=ones(1,nfin-2);
else
    gin=fin(3:nfin)-fin(1:nfin-2);
    gin=2*(gin+(gin==0)).^(-1); % Gain of input triangles
end
msk=fin(2:end-1)==0;
if lowex
    gin(msk)=2*gin(msk);  % double DC input to preserve its power
end
foutin=[fout fin];
nfall=length(foutin);
wleft=[0 fout(2:mfout)-fout(1:mfout-1) 0 fin(2:nfin)-fin(1:nfin-1)]; % left width
wright=[wleft(2:end) 0]; % right width
ffact=[0 gout 0 0 gin(1:min(nf,nfin-nf-2)) zeros(1,max(nfin-2*nf-2,0)) gin(nfin-nf-1:nfin-2) 0]; % gain of triangle posts
% ffact(wleft+wright==0)=0; % disable null width triangles shouldn't need this if all frequencies are distinct
[fall,ifall]=sort(foutin);
jfall=zeros(1,nfall);
infall=1:nfall;
jfall(ifall)=infall; % unsort->sort index
ffact(ifall([1:max(jfall(1),jfall(mfout+1))-2 min(jfall(mfout),jfall(nfall))+2:nfall]))=0;  % zap nodes that are much too small/big

nxto=cumsum(ifall<=mfout);
nxti=cumsum(ifall>mfout);
nxtr=min(nxti+1+mfout,nfall);  % next input node to the right of each value (or nfall if none)
nxtr(ifall>mfout)=1+nxto(ifall>mfout); % next post to the right of opposite type (unsorted indexes)
nxtr=nxtr(jfall);  % next post to the right of opposite type (converted to unsorted indices) or if none: nfall/(mfout+1)
% each triangle is "attached" to the node at its extreme right end
% the general result for integrating the product of two trapesiums with
% heights (a,b) and (c,d) over a width x is (ad+bc+2bd+2ac)*w/6
%
% integrate product of lower triangles
msk0=(ffact>0);
msk=msk0 & (ffact(nxtr)>0); % select appropriate triangle pairs (unsorted indices)
ix1=infall(msk); % unsorted indices of leftmost post of pair
jx1=nxtr(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix1)-foutin(jx1-1); % length of right triangle to the left of the left post
yx=min(wleft(ix1),vfgx); % integration length
wx1=ffact(ix1).*ffact(jx1).*yx.*(wleft(ix1).*vfgx-yx.*(0.5*(wleft(ix1)+vfgx)-yx/3))./(wleft(ix1).*wleft(jx1)+(yx==0));
% integrate product of upper triangles
nxtu=max([nxtr(2:end)-1 0],1);
msk=msk0 & (ffact(nxtu)>0);
ix2=infall(msk); % unsorted indices of leftmost post of pair
jx2=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix2+1)-foutin(jx2); % length of left triangle to the right of the right post
yx=min(wright(ix2),vfgx); % integration length
yx(foutin(jx2+1)<foutin(ix2+1))=0; % zap invalid triangles
wx2=ffact(ix2).*ffact(jx2).*yx.^2.*((0.5*(wright(jx2)-vfgx)+yx/3))./(wright(ix2).*wright(jx2)+(yx==0));
% integrate lower triangle and upper triangle that ends to its right
nxtu=max(nxtr-1,1);
msk=msk0 & (ffact(nxtu)>0);
ix3=infall(msk); % unsorted indices of leftmost post of pair
jx3=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix3)-foutin(jx3); % length of upper triangle to the left of the lower post
yx=min(wleft(ix3),vfgx); % integration length
yx(foutin(jx3+1)<foutin(ix3))=0; % zap invalid triangles
wx3=ffact(ix3).*ffact(jx3).*yx.*(wleft(ix3).*(wright(jx3)-vfgx)+yx.*(0.5*(wleft(ix3)-wright(jx3)+vfgx)-yx/3))./(wleft(ix3).*wright(jx3)+(yx==0));
% integrate upper triangle and lower triangle that starts to its right
nxtu=[nxtr(2:end) 1];
msk=msk0 & (ffact(nxtu)>0);
ix4=infall(msk); % unsorted indices of leftmost post of pair
jx4=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix4+1)-foutin(jx4-1); % length of upper triangle to the left of the lower post
yx=min(wright(ix4),vfgx); % integration length
wx4=ffact(ix4).*ffact(jx4).*yx.^2.*(0.5*vfgx-yx/3)./(wright(ix4).*wleft(jx4)+(yx==0));

% now create the matrix

iox=sort([ix1 ix2 ix3 ix4;jx1 jx2 jx3 jx4]);
msk=iox(2,:)<=(nfall+mfout)/2;
iox(2,msk)=(nfall+mfout+1)-iox(2,msk);  % convert negative frequencies to positive
if highex
    iox(1,iox(1,:)==mfout-1)=mfout-2; % merge highest two output nodes
end
if lowex
    iox(1,iox(1,:)==2)=3; % merge lowest two output nodes
end

x=sparse(iox(1,:)-1-lowex,max(iox(2,:)-nfall+nf+1,1),[wx1 wx2 wx3 wx4],p,nf);

end