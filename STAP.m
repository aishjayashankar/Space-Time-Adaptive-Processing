clc;
clear all;
close all;
count=0;
d(1:199,1)=0;
c = physconst('LightSpeed'); % speed of light;
fc =500; % frequency

lam = c/fc; % wavelength
elementPos = (0:lam/2:lam);
angle = [30;0]; % direction of arrive
sv = steervec(elementPos/lam,angle);
ha = phased.ULA('NumElements',3);
svt=[sv sv sv];
count=0;
Rangebins=199;
NumberofPulses=3;
NumberofAntennas=3;
disp(sv);
Datacube(1:(Rangebins*NumberofPulses),1:(NumberofAntennas)) = rand(Rangebins*NumberofPulses,NumberofAntennas);
Datacube(4,1)=3;
Datacube(203,1)=3;
Datacube(402,1)=3;
Datacube(5,2)=3;
Datacube(204,2)=3;
Datacube(403,2)=3;
Datacube(6,3)=3;
Datacube(205,3)=3;
Datacube(404,3)=3;
Datacube(9,1)=3;
Datacube(208,1)=3;
Datacube(407,1)=3;
Datacube(10,2)=3;
Datacube(209,2)=3;
Datacube(408,2)=3;
Datacube(11,3)=3;
Datacube(210,3)=3;
Datacube(409,3)=3;
disp(Datacube);
RCUT(1:NumberofPulses,1:NumberofAntennas)=0;
y(1:NumberofPulses,1:NumberofAntennas)=0;
x(1:9,1:9)=0;
S(1,1:9)=0;
q(1:9,1:9)=0;
for CellUnderTest=2:196
for i=1:NumberofPulses
RCUT(i,:)=[Datacube(CellUnderTest,i),Datacube(CellUnderTest+Rangebins,i),Datacube(CellUnderTest+(2*Rangebins),i)];
end
for k=1:199
    for j=1:3
    R(j,1)=Datacube(k,j);
    end
    for j=1:3
    R(j,2)=Datacube(k+Rangebins,j);
    end
    for j=1:3
    R(j,3)=Datacube(k+(2*Rangebins),j);
    end
y=R.*svt;
  n=reshape(y,9,1);
  rc=conj(n);
S=transpose(n);
x=rc*S
q=q+x;
end
disp(q);
for k=CellUnderTest-1:1:CellUnderTest+1
for j=1:3
    R(j,1)=Datacube(k,j)
    end
    for j=1:3
    R(j,2)=Datacube(k+Rangebins,j);
    end
    for j=1:3
    R(j,3)=Datacube(k+(2*Rangebins),j);
    end
  y=R.*svt;
  n=reshape(y,9,1);
  rc=conj(n);
S=transpose(n);
x=rc*S;
q=q-x;
end
q=q/196;
%using sv decomposition:
%1. using direct built-in function
%[usvd,ssvd,vsvd]=svd(q);
%Z=vsvd*inv(ssvd)*usvd'

%2. Not using the direct built-in function
q1=q'*q;

[vsvd,ssvd]=eig(q1,'balance');
ssvd1=ssvd.^(1/2);

usvd=q*vsvd;
usvd1=usvd*inv(ssvd1);

svd=usvd1*ssvd1*vsvd';
Z=vsvd*inv(ssvd1)*usvd1'  %<-- Z is the value of inverse of matrix q (verified with original inverse)
%Z=q^-1

fd=(2*1000*5000)/lam;
if  CellUnderTest == 9
    dstv = dopsteeringvec(0,3,500);
else if CellUnderTest==10
    dstv = dopsteeringvec(0,3,500);
 else if CellUnderTest==11
    dstv = dopsteeringvec(0,3,500);
    else
         dstv = dopsteeringvec(fd,3,500);
    end
    end
end
    
B = sv.*dstv;
A=[B;B;B];
C=conj(A);
disp(A);
W=Z*C;
W=transpose(W);
disp(W);
i=reshape(RCUT,9,1);
disp(i);
Thre=W*i;
disp(abs(Thre));
d(CellUnderTest,1)=abs(Thre);
fprintf('Threshold value = %d\n',abs(Thre));
SNRTHRESH = npwgnthresh(1e-6,3,'coherent','linear');
fprintf('NP Threshold value = %d\n',SNRTHRESH);
if(abs(Thre)>(SNRTHRESH))
   fprintf('target is present in rangbin %d',CellUnderTest);
   count =count+1;
else
fprintf('target is not present in rangbin %d',CellUnderTest);
end
end
stem(d);
fft(d);

