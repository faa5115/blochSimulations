

%you want to control the duration, the thickness, excitation profile...
RFDur = 800e-6; %s %original
gamma = 2*pi*42.577*10^6; %rad/s/T. Gyromagnetic ratio.
TBW = 3;%original
%TBW = 4; 
NumOfTPRF = 128; %for the RF profile. 
tRF = linspace(-RFDur/2,RFDur/2,NumOfTPRF); %s.
sliceThickness = 6.00; %mm
flipAngle = 60 * pi/180; %radians. 

%sliceFOV = 30; %mm.
%IntElementsPerMM = 10;
%let's say we want Ns subslices...


Ns = 20;  %Ns must be an integer multiple of sliceThickness.
%so now...our IntElementsPerMM is: 
IntElementsPerMM = Ns/sliceThickness;
subSliceThickness = sliceThickness/Ns; %mm/sub slice

TR = 4e-3;   %s
TE = TR/2; %s

T2 = 150e-3;  %s
%T2 = 1000e-3;
T1 = 1000e-3;  %s

%T1 = 500e-3; %s
%T2 = 500e-3; %s

%T2Fac = 4; %unitless.  just saying after how many T2s will we  begin ignoring signal. 
T2Fac = 3;

%dS is going to be a parameter that defines how quickly the spins flow
%between TRs.  dSMin is always going to be one element shift / TR
%therefore it'll always be 1/Ns.  
%otherwise ... dS will be zero (stationary spins).  
%dS * Ns is the number of elements an isochromat shifted between TRs.  
dSMin = 1/Ns; %smallest nonzero spin exchange percentage.  always going to be 1/Ns.

NPart = 6;
NEffecgtiveFOVElements = Ns * NPart;

%{
dS = 0.3; %fraction. 
%dS = 0;
NShift = round(dS * Ns);
%}
someTempdsValue = 1.0;
NOs =  round(T2Fac  * T2 * someTempdsValue * Ns /TR);
sliceFOV = (Ns + NOs)*subSliceThickness; %mm
%sliceFOV = (NEffecgtiveFOVElements + 2*NOs)*subSliceThickness; %mm


sliceDir = (linspace(-sliceFOV/2, sliceFOV/2, IntElementsPerMM * sliceFOV))'; %mm

%intElementsOnSlice = 100; %number of elements along one axis of a slice.  


RFBW = TBW/RFDur; %Hz %this is the RF's bandwidth. ...
Gslice = RFBW/(gamma*sliceThickness/(2*pi)); %T/mm

B1e = getSincEnvelope_fa(gamma, RFDur, NumOfTPRF, TBW, flipAngle);
%B1efirst = getSincEnvelope_fa(gamma, RFDur, NumOfTPRF, TBW, flipAngle/2);
maxB1e = max(B1e);
%B1e = ones(1, length(B1e));
sliceOffCenter = -4;
%RFPhase0 = pi/2;
RFPhase0 = 0;

%{
NumOfExc = 300;
M0     = zeros(size(sliceDir,1), 3);
MbSSFP = zeros(size(sliceDir,1), 3, NumOfExc );

Meq = zeros(size(sliceDir,1), 1);
for n = 1 : size(sliceDir, 1)
    M0(n, :) = [0 0 1]';
    Meq(n)   = norm(squeeze(M0(n,:) ));
end
%}
phaseInc = pi;

NumOfOffRes = 75;
PhiPerTR = linspace(-2*pi, 2*pi, NumOfOffRes);
OffRes   = PhiPerTR/(2 * pi * TR);

offresonanceones = ones(size(sliceDir));

%% Now we will include pulsatile flow. 
%lets first define one flow period. 
%let's make it a multiple of our TR...to make our life easier. 
%our bloch simulation earlier had a TR = 4 ms. 
%let's make a period be 900 ms. 

numOfTP_og = 250;
%subTractFromTP_og = mod(numOfTP_og, intendedVPS);
numOfTP = numOfTP_og;% - subTractFromTP_og;

NumOfExcPerBeat = numOfTP;
t_pulse = 0 : TR : (NumOfExcPerBeat - 1) * TR; %ms
vMax = 1500; %mm/s
dSMax = TR * vMax/sliceThickness; 
%PeakVelExc = round(NumOfExcPerBeat / 3); 
PeakVelExc = round( numOfTP * 0.18 );
stdDevExc = 10;  
dSsim_pulse = dSMax * exp(- ((t_pulse - PeakVelExc * TR).^2)/(2*((stdDevExc*TR)^2)));
NumOfHeartBeats = 4;

dSsim_total = zeros(1, NumOfHeartBeats * length(dSsim_pulse) );

for n = 1 : NumOfHeartBeats
    dSsim_total(1, (n - 1) * NumOfExcPerBeat + 1 : n * NumOfExcPerBeat  ) = dSsim_pulse;
end



NumOfExc = NumOfExcPerBeat * NumOfHeartBeats;
M0     = zeros(size(sliceDir,1), 3);
MbSSFP = zeros(size(sliceDir,1), 3, NumOfExc );

Meq = zeros(size(sliceDir,1), 1);
for n = 1 : size(sliceDir, 1)
    M0(n, :) = [0 0 1]';
    Meq(n)   = norm(squeeze(M0(n,:) ));
end

excitationArr = 1 : length(dSsim_total);
%%
OffResIter = round(length(PhiPerTR)/2);
        df = squeeze(OffRes(1,OffResIter));
        [Atr, Btr] = freeprecess(TR,T1,T2,df);

        MSliceProfile0 = func_sliceSelection (B1e, sliceDir, RFBW, sliceThickness, RFDur, M0, gamma, sliceOffCenter, RFPhase0  , Meq, T1, T2, df *  offresonanceones );

        MSliceProfile=func_sliceSelRef (RFDur, RFBW, sliceThickness,  B1e, squeeze(MSliceProfile0(:, :, size(MSliceProfile0,3))), sliceDir, Meq, T1, T2, df *  offresonanceones);
        %{
        figure, plot( sliceDir,angle (squeeze(MSliceProfile(:,1,size(MSliceProfile,3))) +1i* squeeze(MSliceProfile(:,2,size(MSliceProfile,3)))))

        figure, plot( sliceDir,abs ( squeeze(MSliceProfile(:,3,size(MSliceProfile,3)))))


        figure, cplot(  (squeeze(MSliceProfile(:,1,size(MSliceProfile,3))) +1i* squeeze(MSliceProfile(:,2,size(MSliceProfile,3))))  )
        %}

        % calculate the flip angle and the phase of each isochromat.   --------------------------------------------

        localNetRFPhase = zeros(length(sliceDir), 1);
        localNetRFFA    = zeros(length(sliceDir), 1);

        %in the loop below ... you changed the signs and swapped the im/re
        %components because rf phase is 90 deg out of phase of the magnetization's
        %phase.  
        for n = 1: length(sliceDir)
            %fadilali test.  gives the proper rephasing at TE = TR/2.
            localNetRFPhase(n, 1) = angle(  squeeze(MSliceProfile(n,1,size(MSliceProfile,3))) + 1i *squeeze(MSliceProfile(n,2,size(MSliceProfile,3))));
            %localNetRFPhase(n, 1) = -pi/2;
            localNetRFFA   (n, 1) = angle( 1i*abs(squeeze(MSliceProfile(n,1,size(MSliceProfile,3))) +1i* squeeze(MSliceProfile(n,2,size(MSliceProfile,3)))) +  (squeeze(MSliceProfile(n,3,size(MSliceProfile,3)))  )     );
        end
        % calculate the proper rotation matrix for each isochromat.  
        estExcMat = zeros(3, 3, length(sliceDir)); %estimated excitation matrix. 
        for n = 1 : length(sliceDir)
            estExcMat(:, :, n) = throt(squeeze ( localNetRFFA(n, 1) ), squeeze( localNetRFPhase(n, 1) ) );
        end

        localNetRFPhase0 = localNetRFPhase;
        localNetRFFA0 = localNetRFFA;
        estExcMat0 = estExcMat;
MAllExcitationsEachOffResCycle = zeros(size(M0,1), size(M0,2), NumOfExc, NumOfOffRes);


%%

MEndExc = squeeze(MSliceProfile(:, :, squeeze(size(MSliceProfile,3)) ));
figure, 
hold on

yyaxis left
plot(sliceDir, abs(squeeze(MEndExc(:, 1)) + 1i*squeeze(MEndExc(:, 2))), 'linewidth', 5.0 )
yyaxis right
 ylim([-pi, pi])
plot(sliceDir, angle(squeeze(MEndExc(:, 1)) + 1i*squeeze(MEndExc(:, 2))), 'linewidth', 5.0 )
legend('Magnitude', 'Phase')
xlim([-10 10])
title('Transverse Magnetization:  Magnitude and Phase')
set(gca,'FontSize',20, 'FontWeight', 'Bold')
xlabel('slice location (mm)')

%% plot b1e profile
rftimearr = linspace(0, RFDur, NumOfTPRF);
B1ea = zeros(size(B1e));
tRF = linspace(-RFDur/2,RFDur/2,NumOfTPRF); %s.
for n = 1 : size(B1e,2)
    B1ea(1,n) = B1e(1,n) * (cos(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n))) + 1i*(sin(RFPhase0 + gamma*Gslice*sliceOffCenter*(tRF(1,n)))));
end

figure,
cplot2(rftimearr, B1ea )
legend('Real', 'Imaginary')
set(gca,'FontSize',20, 'FontWeight', 'Bold')
xlabel('time (ms)',  'fontsize', 20 )
ylabel('B1 Magnitude')
