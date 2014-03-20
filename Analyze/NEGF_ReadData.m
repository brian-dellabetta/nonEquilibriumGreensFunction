clear all; format long;

%insert path to output folder here
Folder = '..\TRAB1Corr2D0Dis3D-1_Mat0_a01D0_30x10x6_bc000_V-1D0_0D0_muT0D0_ABPhi5D-1';
cd (Folder);

q = 1.602D-19; h = 6.626D-34; kB = 1.38066D-23; a1 = 2.456D-10; b0 = 3.35D-10; T0=0.0; Vt=kB*T0/q;
eps0 = 8.85418D-12; m0 = 9.1019D-31; hbar = h/(2D0*pi); 
mu0 = pi*4e-7; alpha = q^2/hbar; c0 = 1;
labelFont = 20; axisFont = 16;
%Natural Units m0=1.0;meff=1.0; q=1.0;hbar=2.0*pi;hbar=1.0;kB=1.0;

results = load('scalars.txt'); cursor = 1;
Npx = results(cursor); cursor = cursor+1;
Npy = results(cursor); cursor = cursor+1;
Npz = results(cursor); cursor = cursor+1;
Norb = results(cursor); cursor = cursor+1;
muT = results(cursor); cursor = cursor+1;
eMin = results(cursor); cursor = cursor+1;
eMax = results(cursor); cursor = cursor+1;
eStpnum = results(cursor); cursor = cursor+1;
deltaE = results(cursor); cursor = cursor+1;
mu = zeros(1, 2);
for ii = 1:2
    mu(ii) = results(cursor);
    cursor = cursor+1;
end
channelType = results(cursor); cursor = cursor+1;
if (channelType == 1)
    a0 = results(cursor); cursor = cursor+1;
else
    a0 = 1.0;
end
abPhi = results(cursor); cursor = cursor+1;
corrLength = results(cursor); cursor = cursor+1;
DisStr = results(cursor); cursor = cursor+1;

NN = Norb*Npy; TNp = Norb*Npx*Npy*Npz;   %Total number of points
NNx = Norb*Npy*Npx; NNz = Norb*Npy*Npz;

%allocate arrays to be read in
RhoN = zeros(Npx,Npy,Npz);
RhoP = zeros(Npx,Npy,Npz);

EDensity = zeros(eStpnum,1);
HDensity = zeros(eStpnum,1);
ContactCurrent = zeros(1,2);
Ip = zeros(2,1);
T12 = zeros(eStpnum, 1);
Modec1 = zeros(eStpnum, 1);
Modec2 = zeros(eStpnum, 1);
EAxis = zeros(eStpnum, 1);
DOSc1 = zeros(eStpnum, 1);
DOSc2 = zeros(eStpnum, 1);
Jx = zeros(Npx, Npy, Npz);
Jy = zeros(Npx, Npy, Npz);
Jz = zeros(Npx, Npy, Npz);

results = load('Greens.txt'); cursor = 1;
for xx = 1:Npx
    for yy = 1:Npy
        for zz = 1:Npz
            RhoN(xx,yy,zz) = results(cursor); cursor = cursor+1;
            RhoP(xx,yy,zz) = results(cursor); cursor = cursor+1;
        end
    end
end
if (cursor+Npx < length(results))
    for xx = 1:Npx
        for yy = 1:Npy
            for zz = 1:Npz
                DeltaR(xx,yy,zz) = results(cursor);
                cursor = cursor+1;
            end
        end
    end
end

results = load('Transmissions.txt'); cursor = 1;
for ii = 1:eStpnum
    EDensity(ii) = results(cursor); cursor = cursor+1; Npx*Npy*Npz;
    HDensity(ii) = results(cursor); cursor = cursor+1; Npx*Npy*Npz;
end
for ii = 1:eStpnum
    T12(ii) = results(cursor); cursor = cursor+1;
end
for ii = 1:eStpnum
    Modec1(ii) = results(cursor); cursor = cursor+1;
    Modec2(ii) = results(cursor); cursor = cursor+1;
end
for ii = 1:eStpnum
    EAxis(ii) = results(cursor);
    cursor = cursor+1;
end
for ii = 1:2
    ContactCurrent(ii) = results(cursor)*q;
    cursor = cursor+1;
end
for xx = 1:Npx
    for yy = 1:Npy
        for zz = 1:Npz
            Jx(xx, yy, zz) = results(cursor);
            cursor = cursor+1;
            Jy(xx, yy, zz) = results(cursor);
            cursor = cursor+1;
            Jz(xx, yy, zz) = results(cursor);
            cursor = cursor+1;
        end
    end
end
for ii = 1:eStpnum
    DOSc1(ii) = results(cursor);
    cursor = cursor+1;
    DOSc2(ii) = results(cursor);
    cursor = cursor+1;
end

Jtot = sqrt((Jx).^2+(Jy).^2+(Jz).^2);
LDOS = EDensity+HDensity;

phi = zeros(Npx,Npy,Npz);
ro = zeros(Npx,Npy,Npz);
results = load('Poissons.txt'); cursor = 1;
if (length(results) > Npx*Npy*Npz)
    for xx = 1:Npx
        for yy = 1:Npy
            for zz = 1:Npz
                phi(xx, yy, zz) = results(cursor);
                cursor = cursor+1;
                ro(xx, yy, zz) = results(cursor);
                cursor = cursor+1;
            end
        end
    end
end

if (exist('SigExp.txt'))
    results = load('SigExp.txt'); cc = 1;
    if (length(results) > Npx*Npz)
        SigXExp = zeros(Npx,Npy,Npz);
        SigYExp = zeros(Npx,Npy,Npz);
        SigZExp = zeros(Npx,Npy,Npz);
        for xx = 1:Npx
            for yy = 1:Npy
                for zz = 1:Npz
                    SigXExp(xx,yy,zz) = results(cc); cc = cc+1;
                    SigYExp(xx,yy,zz) = results(cc); cc = cc+1;
                    SigZExp(xx,yy,zz) = results(cc); cc = cc+1;
                end
            end
        end
    end
end
    
results = load('ERCD.txt'); cursor = 1;
if(length(results)>Npy*Npx*Npz*eStpnum)
    ERJx = zeros(Npx,Npy,Npz,eStpnum);
    ERJy = zeros(Npx,Npy,Npz,eStpnum);
    ERJz = zeros(Npx,Npy,Npz,eStpnum);
    ERRhon = zeros(Npx,Npy,Npz,eStpnum);
    for xx = 1:Npx
        for yy = 1:Npy
            for zz = 1:Npz
                for ee = 1:eStpnum
                    ERJx(xx, yy, zz, ee) = results(cursor);
                    cursor = cursor+1;
                    ERJy(xx, yy, zz, ee) = results(cursor);
                    cursor = cursor+1;
                    ERJz(xx, yy, zz, ee) = results(cursor);
                    cursor = cursor+1;
                    ERRhon(xx, yy, zz, ee) = results(cursor);
                    cursor = cursor+1;
                end
            end
        end
    end
end
clear results;

figure(1); hold on; plot(EAxis, T12, 'k'); hold on;
set(gcf, 'Position', [10 500 570 360]); %set(gca, 'FontSize', axisFont, 'FontWeight', 'bold'); 
title('Transmission (e^2/h)', 'fontsize', labelFont, 'FontWeight', 'bold'); 
xlabel('Energy (eV)', 'fontsize', labelFont, 'FontWeight', 'bold'); 
ylabel('Transmission (e^2/h)', 'fontsize', labelFont, 'FontWeight', 'bold'); 


figure(2); set(gcf, 'Position', [560 500 570 360]);
quiver(squeeze(Jx(:,:,1))', squeeze(Jy(:,:,1))'); 

cd ..


