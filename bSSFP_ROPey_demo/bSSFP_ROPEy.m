P = phantom('Modified Shepp-Logan',80);
%%
Gmax = 40e-3; %T/mm 

T1 = 1170e-3; %s
T2 = 66e-3;   %s
%T1 = 500e-3; %s
%T2 = 500e-3; %s

NPEy = 60;
NRO  = 60; 

lengthAlongPEy = 50; %mm.
lengthAlongRO  = 50; %mm.
%The actual object, P, has physical dimensions lengthAlongPEy x lengthAlongRO.

%IntElementsPerMMy = round(size(P,2)/lengthAlongPEy);
%IntElementsPerMMx = round(size(P,1)/lengthAlongRO);

IntElementsPerMMy = (size(P,2)/lengthAlongPEy);
IntElementsPerMMx = (size(P,1)/lengthAlongRO);

lengthPEy = (linspace(-lengthAlongPEy/2, lengthAlongPEy/2, IntElementsPerMMy * lengthAlongPEy))'; %mm
lengthRO = (linspace(-lengthAlongRO/2, lengthAlongRO/2, IntElementsPerMMx * lengthAlongRO))'; %mm

Object = zeros(size(P,1), size(P,2), 3);
Object(:, :, 3) = P;
%KSpace = zeros(NRO, NPEy, 3);

gamma = 2*pi*42.577*10^6; %rad/s/T. Gyromagnetic ratio.

BWpp = 500; %Hz/pixel.  ReadOut bandwidth per pixel. 
dwellTime = 1/(BWpp * NRO);
FOVro  = 50; %mm
FOVpey = 50; %mm
dKro = 1/FOVro;
dKpey = 1/FOVpey;


Gro = dKro * 2 * pi * (BWpp * NRO) / gamma; %Readout gradient. T/mm. 

PhaseROPPMM = gamma * round(NRO/2) * dwellTime * -Gro; %phase/mm  caused by
                                                       %the readout 
                                                       %prephaser
%unit check:
%(rad/T/s) * (s) * (T/mm) = rad/mm. 

%calculate the max time needed for the phase encoding gradient...
%assuming we ramp up to the max gradient amplitude ...
%you know what...fuck it.  because i'm lazy, i'll just make the PEy and 
%ROP have the same duration. 

%let's calculate the minimum TE.  we'll use that TE.  and of course... our 
%TR will be 2* TE. 
%oh...and we already have our 2D slice (the shepp-logan phantom) ... so
%we'll just use standard matrix rotation for the excitation. 
%we'll use the dwellTime as our increments.  because i'm too lazy to reraster. 

%first time point is for the excitation. 
%TEmin = dwellTime;
%then PEy and ROP:
%TEmin = TEmin + round(NRO/2) * dwellTime;
TEmin = round(NRO/2) * dwellTime;
%then the readout's first half:
TEmin = TEmin + round(NRO/2) * dwellTime;

%sooooo....
TE = TEmin;
TR = 2*TE; 



prepTRs = 500;
MPrep = zeros(size(Object,1), size(Object,2), size(Object,3), prepTRs);
%MPrep(:, :, :, 1) = Object(:, :, :);

df = 0;
%let's have off resonance varying by space an option. 
%ATR = zeros(size(P,1),size(P,2), 3, 3);
%BTR = zeros(size(P,1),size(P,2), 3, 1);

%{
for n = 1 : size(P,1)
    for m = 1: size(P,2)
        %off resonance slope factors along x or y. .. make them zero for now.
        dFx = 0;
        dFy = 0;
        
        df = dFx * (lengthRO(n,1)) + dFy * (lengthPEy(m,1));
        
        [ATR(n, m, :, :), BTR(n, m, :, :)] = freeprecess(TR,T1,T2,df);
    end
end
%}

%off resonance slope factors along x or y. .. make them zero for now.
dFx = 0;
dFy = 0;
        
flipAngle = pi/2;%flip angle. radians.  
phase = 0;  %RF phase. radians. 
phaseInc = pi;

%just to have each spin reach steady state. 
for k = 0 : prepTRs - 1
    exciteMatrix = throt(flipAngle, phase);
    phase = phase + phaseInc;
    for n = 1 : size(P,1)
        for m = 1 : size(P,2)
            
            %MPrep(n, m, :, p ) holds the magnetization after p excitations
            %for the location (n, m).  
            
            if (k+1 == 1)
                MPrep(n, m, :, k+1) = exciteMatrix * squeeze(Object(n, m, :)) ;
            else
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                dfGradDeph = 2 * pi * TR * lengthRO ;
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                df = dFx * (lengthRO(n,1)) + dFy * (lengthPEy(m,1)) + squeeze(dfGradDeph(n, 1));
                [ATR, BTR] = freeprecess(TR,T1,T2,df);
                MPrep(n, m, :, k+1) = exciteMatrix * (ATR * squeeze(MPrep(n, m, :, k)) + BTR * norm(squeeze(Object(n, m, :))) );
            end
            
            % so...
            % MPrep(x1, y1, :, k) holds the magnetization immediately after k
            % excitations at location (x1, y1). It does NOT pursue the decay afterwards. 
        end
    end
    disp(k)
end
%%
MPrepxy = squeeze(MPrep(:, :, 1, : ) + 1i * MPrep(:, :, 2, : ));

%%
figure, plot(abs( squeeze(MPrepxy(44, 72, :)) ) ) 

%%
figure,
imagesc(abs(squeeze(MPrepxy(:, :, size(MPrepxy,3) )))),
colormap('gray')
%%
figure,
imagesc(abs(Object(:, :, 3))),
colormap('gray')
%% now that we've reached steady state... let's do an imaging readout. 
%exciteMatrix * (ATR * squeeze(MPrep(n, m, :, k)) + BTR * norm(squeeze(Object(n, m, :))) );
%MStart = zeros(size(MPrep, 1 ) , size(MPrep, 2), size(MPrep, 3));

%just to have the full TR decay.
%{
for n = 1 : size(MPrep, 1)
    for m = 1 : size(MPrep, 2)
        df = dFx * (lengthRO(n,1)) + dFy * (lengthPEy(m,1));
        [ATR, BTR] = freeprecess(TR,T1,T2,df);
        MStart(n, m, :) = exciteMatrix * (ATR * squeeze(MPrep(:, :, :, size(MPrep,4))) + BTR * norm(squeeze(Object(n, m, :))) );
    end
end
%}
MStart = squeeze(MPrep(:, :, :, size(MPrep,4)));
%---------------------------------------FADIL ALI----------------------
% When you return... you must debate to see if you want to start with an
% excitation in the loop below... or if you want to end with the
% excitation...if you want to end with an excitation ... then you should
% return to your original MStart above...
%-----------------------------------------------------------------------


MPreph = zeros(size(MStart));
angleROP = PhaseROPPMM * lengthRO;
%calculate the rotatoin about z due to the RO prephaser.
zrotROP = zeros(3, 3, length(lengthRO));
dGpey = dKpey * 2 * pi / (gamma * round(NRO/2) * dwellTime ) ;
for n = 1 : length(lengthRO)
    zrotROP(:, :, n) = zrot(squeeze(angleROP(n)));
end


zrotPEy = zeros(3, 3, length(lengthPEy));
peInt = 1; 
MReadOutData = zeros(size(Object,1), size(Object,2), 3, NRO, NPEy);
%remember...as we enter this loop...we cannot lose the steady state. 


phaseArr = zeros(NPEy, 1);
for pestep = round(NPEy/2) : -1 : -1 * round(NPEy/2) + 1  %iterate through each phase encoding step
    %disp(peInt)
    disp(pestep)
    phaseArr(peInt,1) = phase; 
    exciteMatrix = throt(flipAngle, phase);
    %just to make sure there is no added phase offsent between excitations
    %from the RF's phase. 
    %rotRFPhaseConsideration = zrot(-phase); 
    %rotRFPhaseConsideration = zrot(phase); 
    phase = phase + phaseInc;
    %first calculate the phase encoding phase disperion.  
    Gpey = pestep * dGpey; 
    disp(Gpey)
    
    PhasePEyPMM = gamma * Gpey * round(NRO/2) * dwellTime ; %rad/mm. 
    %untie test:  rad/T/s * T/mm * s = rad/mm. 
    %above is accurate:
    %(rad/T/s) * (T/mm) * (s) = rad/mm.
    anglePEy = squeeze(PhasePEyPMM) * lengthPEy; 
  
    %there does not seem to be a problem with this part...
    %figure, plot(anglePEy), title(strcat('pe ', num2str(pestep), 'PhasePEyPMM ', num2str(PhasePEyPMM), ' Gpey ', num2str(Gpey)  ))
    
    %we are going to iterate through all space and change the magnetization
    %vector at each location according to the ROP, PEy, df phase
    %dispersions...and through T2 decay and T1 recovery. 
    
    %FADILALI CHECKPOINT HERE FOR REFERRING TO THIS FILE FOR CONSTRUCTING
    %YOUR VOLUMETRIC MOTION CODE!!!!
    
    for n = 1 : size(MStart, 1) 
        for m = 1 : size(MStart, 2)
            df = squeeze(dFx * (squeeze(lengthRO(n,1))) + dFy * (squeeze(lengthPEy(m,1))));
            [APreph, BPreph] = freeprecess(round(NRO/2) * dwellTime, T1,T2,df);
             zrotPEy(:, :, m ) = zrot(squeeze(anglePEy(m)));
             
             %calculate the magnetization vector at each location as a
             %result of the PEy gradient and the Readout prephaser. 
                                                       
 MPreph(n, m, : ) = APreph * squeeze(zrotPEy(:, :, m))*squeeze(zrotROP(:, :, n)) * squeeze(MStart(n, m, :)) + BPreph * norm(squeeze(Object(n, m, :)));
        end
    end
    %debug
    %this part seems to be good...
    %figure, imagesc(angle( squeeze(MPreph(:, :, 1)) + 1i*squeeze(MPreph(:, :, 2)) ) ) ,colormap('gray'), title(strcat('pe ', num2str(pestep), 'PhasePEyPMM ', num2str(PhasePEyPMM), ' Gpey ', num2str(Gpey)  ))
    
    %soo..above we iterated through each isochromat to prephase and phase
    %encode them.  
 %okay... so we now just pre encoded everything we need...now for the
 %readout.  we are going to be sampling at each dwell time...so we'll need
 %the appropriate decay matrix for each...at eac hlocation (because of off
 %resonance varying spatially).  of the format: 
 %           [Adt, Bdt] = freeprecess(dwellTime, T1,T2,df);
 
 %for each timepoint within the readout... we must go through each position
 %in space to calculate the magnetization rotation.  each rotation will
 %depend on the previous timepoint. 
    roInt = 1;

    %debugging when commented out.
    %so now...we will iterate through each isochromat ... to record how its
    %magnetization precesses at each time point of the readout.  
    for n = 1 : size(MStart, 1) 
        for m = 1 : size(MStart, 2)
            for roInt = 1 : NRO
                if (roInt == 1)
                    MBefore =  squeeze(MPreph(n,m, :));
                else
                    MBefore = squeeze(MReadOutData(n, m, :, roInt - 1, peInt));
                end
                
                df = dFx * squeeze(lengthRO(n,1)) + dFy * squeeze(lengthPEy(m,1));
                df = df + gamma * Gro * squeeze(lengthRO(n))/(2 * pi);
                %unit check:  (rad/T/s) * (T/mm) * (mm) * (1/rad) = 1/s.  
                [Adt, Bdt] = freeprecess(dwellTime, T1,T2,df);
                
                MReadOutData(n, m, :, roInt, peInt ) = Adt  * squeeze(MBefore) + Bdt * norm(squeeze(Object(n, m, :)));
            end
        end
    end
    
        
    %by this point all of our spins have gone through the readout... it is
    %now time to balance the readout and the PE lobes.  
     %MTREnd = zeros(size(Object, 1), size(Object, 2), size(Object, 3));
     
     %debugging when commented out
     
     for n = 1 : size(MStart, 1) 
        for m = 1 : size(MStart, 2)
            df = dFx * (lengthRO(n,1)) + dFy * (lengthPEy(m,1));
            %we can use the same prephaser from above.  because it's the
            %same time duration and inherant off resonance. 
            [APreph, BPreph] = freeprecess(round(NRO/2) * dwellTime, T1,T2,df);
            
             PhaseRePEyPMM = gamma * -Gpey * round(NRO/2) * dwellTime ; %rad/mm. 
             angleRePEy = PhaseRePEyPMM * lengthPEy; 
            %phase encode rewinder.  
            
             zrotPEy(:, :, m ) = zrot(squeeze(angleRePEy(m)));
             
             %debug
             %zrotPEy(:, :, m ) = zrot(0);
             
             %use the same rotation matrix for the readout prephaser as
             %above. 
             
             %calculate the magnetization vector at each location as a
             %result of the PEy gradient and the Readout prephaser. 
             %then we also excite after the remaining decay.  
 MStart(n, m, :) = exciteMatrix * (APreph * squeeze(zrotPEy(:, :, m)) * squeeze(zrotROP(:, :, n)) * squeeze(MReadOutData(n, m, :, NRO, peInt)) + BPreph * norm(squeeze(Object(n, m, :))) );
        end
    end
       
    
    peInt = peInt + 1;
end

%%
MPrephxy = squeeze(MPreph(:, :, 1 ) + 1i * MPreph(:, :, 2 ));
%%

figure,
imagesc(angle(squeeze(MPrephxy(:, :)))),
colormap('gray')

%%
n = 30;
m = 30;
figure,
imagesc( angle(squeeze(MReadOutData(:, :, 1, n, m)) + 1i*squeeze(MReadOutData(:, :, 2, n, m))  )),
colormap('gray')
%%
% MReadOutData is 
%  Object's length along readout, object's length along pey, 
% [Mx, My, Mz]', each duration (segmented by dwellTime) along RO, for each
% PE step. 
%KSpace is a vector that holds the complex points during each readout. 
KSpace = zeros(NRO, NPEy);
for roInt = 1 : NRO
    for peInt = 1 : NPEy
        
        %sum up all of the isochromates for that PE step for that
        %particular point in the readout. 
        
        rotRFPhaseConsideration = zrot(-phaseArr(peInt)); 
        
        for n = 1 : size(Object, 1)
            for m = 1 : size(Object, 2)
                MTemp = rotRFPhaseConsideration * squeeze(MReadOutData(n, m, :, roInt, peInt)); 
                KSpace(roInt, peInt) = KSpace(roInt, peInt) + (squeeze(MTemp(1)) + 1i * squeeze(MTemp(2)));
            end
        end
        
    end
    disp(roInt)
end
%%
figure,
imagesc(abs(KSpace)),
colormap('gray')

%%
figure,
imagesc(abs(ifftnc(KSpace))),
colormap('gray')
%%
figure,
imshow(abs(ifftnc(KSpace)), []),

%% k space half.
startPF = 20;
KSpaceHalf = KSpace(:, startPF:end);
KSpaceHalf0Pad = zeros(size(KSpace));
KSpaceHalf0Pad(:, startPF:end) = KSpaceHalf;
%%
figure,
imagesc(abs(fftnc(KSpaceHalf)))
title('nothing')
colormap('gray')
%this just gives poor spatial resolution.  the ky max is shorter. 

figure,
imagesc(abs(fftnc(KSpaceHalf0Pad)))
title('zero pad')
colormap('gray')

figure,
imagesc(abs(fftnc(KSpace))),
colormap('gray')
title('full acquisition')

% notice the Gibbs ringing along the partial fourier direction...