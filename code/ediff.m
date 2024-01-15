% --------------------------------------------------------------------------------
%  Electron diffraction - P. Ahrenkiel
% --------------------------------------------------------------------------------
tic
clear all

%% Flags
doPlotExpt=true;%
doPlotSim=true;%

%% Composition and etaMax
% (Al(w)Ga(1-w))(x)In(1-x)P
xIn=                          0.63;%<-Input

w=1;x=1-xIn;y=0;z=0;
a=a_nm().get_a(w,x,y,z);

% Find eta max
xMin=xIn;
if 1-xIn<xMin
    xMin=1-xIn;
end
etaMax=2*xMin;

%Default values
eta=                        0.3618;   %default order parameter
thick=                      100;    %default thickness

%% Experiments <- Input
expt1=                              Expt('Experiment 1');
%expt1.filename=                     "E:\Documents\MATLAB\ediff\data\3_7191_1PED_A_5_ext.txt";%3-7192-1-A0.txt";%3-7192-1-A0_1.txt";%3-7191-2-4-13_2.txt";%";%data file             
%expt1.filename=                     "E:\Documents\MATLAB\ediff\data\3_7191_1PED_A_5_ext.txt";%3-7192-1-A0.txt";%3-7192-1-A0_1.txt";%3-7191-2-4-13_2.txt";%";%data file             
expt1.filename=                     "E:\Documents\MATLAB\ediff\data\3-7192-1-A0.txt";             
%expt1.filename=                     "E:\Documents\MATLAB\ediff\data\3-7191-2-4-13_2.txt";             
%expt1.filename=                     "E:\Documents\MATLAB\ediff\data\3_7191_1PED_A_5_ext.txt";%3-7192-1-A0.txt";%3-7192-1-A0_1.txt";%3-7191-2-4-13_2.txt";%";%data file             
%expt1.filename=                     "E:\Documents\MATLAB\ediff\data\3_7191_1PED_A_5_ext.txt";%3-7192-1-A0.txt";%3-7192-1-A0_1.txt";%3-7191-2-4-13_2.txt";%";%data file             
expt1.doPlot=                       doPlotExpt;%
expt1.selectData=                   false;
expt1.balanceIntens=                false;
expt1.includePairDiffError=         true;
%{
expt2=                      Expt('Experiment 2');
expt2.filename=             "E:\Documents\MATLAB\ediff\data\3-7191-A5_1.txt";%data file             
expt2.doPlot=               doPlotExpt;%
expt2.selectData=           false;
%}
%% constants
const=                      Const();

% Basis
basis=                      Basis();
basis.A_refM=              	a*[ 1 0 0;
                                0 1 0;
                                0 0 1];%Vectors are rows

%Get inverse basis
basis.B_refM=inv(basis.A_refM)';

%% Conditions <- Input

%Diffraction conditions---------------------------
diff1=                      Diff(const,200);
diff1.T_K=                  300;%Not currently used.

diff1.phi=                  20;%beam tilt (mrad)
diff1.nSteps=               90;%# of precession steps
diff1.N_ZoneL=              [0];%zone index list <- Input
diff1.sgRange=              0.2;%Buerger precession filter range (1/nm)
diff1.mode=                 'max';%maximum: 'max',integrated: 'int'
diff1.doBuerger=            false;%Use Buerger mask for precession

diff1.gMax=                	20/a;%(1/nm) largest g to plot

diff1.tExposure=         	0;%
diff1.fContrast=            15;%
diff1.camlen=           	24000;%pix
diff1.pixX=                 2048;%image size X
diff1.pixY=                 2048;%image size Y

%Get incident wavevector
diff1.kV=                   const.nZ/diff1.lambda;
diff1.doRealPlot=             true;%more realistic plot of data for export

%% Define phases <- Input
%zincblende matrix---------------------------
matrix=                     DisordPhase('Matrix').AlInP(xIn);

%AlInP (-1 1 1) variant---------------------------
varL=                       OrdPhase('AlInP variant L',etaMax).AlInP_1n11(xIn,eta);

%AlInP (1 -1 1) variant---------------------------
varR=                       OrdPhase('AlInP variant R',etaMax).AlInP_n111(xIn,eta);

%Si---------------------------
%Si=                         DisordPhase('Si').Si();

%% Domains <- Input
matrix1=                    Domain(matrix);
ordL1=                  	Domain(varL);
ordR1=                  	Domain(varR);
%wafer=                  	Domain(Si);

%% Regions <- Input
%Region 1---------------------------
region1=                    Region("Region 1",basis);

region1.thick=              thick;%nm%<-Input
region1.coherThickRat=      0.1;%coherent fraction
region1.incoherThickRat=    0.5;%fractional variation in thickness
region1.uvwZoneV=           [1 1 0]';%zone axis%<-Input
region1.hklXV=              [-2 2 0]';%X-drection, specified by RLV<-Input

%Add domains
%region1.addDomain(matrix1,0.4);
region1.addDomain(ordR1,0.5);%.addDomain(ordL1,0.5);%.addDomain(matrix1,0.3);

%% Simulations <- Input
sim1=                       Sim('Simulation 1');   
sim1.gMax=                  16/a;%(1/nm) largest g to consider
sim1.sgMax=                 0.01;%(1/nm) always include beam in diagonalization if |sg|<sgMax
sim1.invwMin=               0.002;%include as perturbation if |sg|>sgMax and 1/|sg*Xsig|>invwMin
sim1.useAllFitBeams=        false;%include all beams satisfying fit conditions, regardless of above.
sim1.usePert=              	false;
sim1.useCoher=              true;
sim1.useIncoher=            true;
sim1.doPlot=                doPlotSim;%<- Input
sim1.diff=                  diff1;
sim1.region=                region1;
sim1.dynamical=             true;
sim1.useReflection=         true;

%% Fitting information <- Input
fit1=                       Fit(expt1,sim1,region1,diff1);
fit1.gMax=                  8.1/a;%(1/nm) g<=gMax to fit
fit1.gMin=                  0./a;%(1/nm) g>=gMin to fit
fit1.plotFitBeamsOnly=      false;%whether to only plot fitted beams<- Input
fit1.fitPairsOnly=          true;
fit1.fitF=                  false;%fit structure factors: true, intensities: false
fit1.H_specRefL=            [%Specific beams to include
%{
                                [	-2	2	-4	]
                                [	-1.5	1.5	4.5	]
%}
                            ];

proj=Proj(const,basis);
proj.addFit(fit1);%.addFit(fit2);

%(Un)comment items from list below to specify which parameters to vary.
proj.parL=                  {
                            {'thickness';'Region 1'}
                            %{'thickness';'Region 2'}
                            {'eta';'AlInP variant R'}
                            %{'eta';'AlInP variant L'}
                            %{'volume';'Region 1';'Matrix'}
                            {'volume';'Region 1';'AlInP variant R'}
                            %{'volume';'Region 2';'AlInP variant R'}
                            };

%(Un)comment items from list below to specify which parameters to vary.
proj.constrL=                {
                            {'eta equal';'AlInP variant L';'AlInP variant R'}
                            %{'eta/volume ratio constant'}
                            %{'eta max';'Variant 1'}
                            %{'eta max';'Variant 2'}
                    	};

%(Un)comment items from list below to specify which variables to report as text.
proj.varL=                  {
                            {'thickness';'Region 1'}
                            {'eta';'AlInP variant R'}
                            %{'eta';'AlInP variant L'}
                            %{'volume fraction';'Matrix'}
                            {'volume fraction';'Region 1';'AlInP variant R'}
                            %{'volume fraction';'Region 1';'AlInP variant L'}
                            {'eta net'}
                            };

%Used for sweeps
proj.sweepL=                {
                            {'thickness';'Region 1';10.^(0:3/40:3)'}
                          	{'eta';'AlInP variant R';(0.0:etaMax/25:etaMax)'}
                           	%{'volume';'Region 1';'AlInP variant R';10.^(-2:3/30:2)'}
                         	};
                    

proj.sweepParL=             {
                            %{'thickness';'Region 1'}
                            %{'eta';'AlInP variant R'}
                            %{'eta';'AlInP variant L'}
                            %{'volume';'Region 1';'Matrix'}
                            %{'volume';'Region 1';'AlInP variant R'}
                            %{'volume';'Region 2';'AlInP variant R'}
                            };

%max. # of iterations
proj.nIters=                 10;
proj.report=                true;
proj.init();

%sim1.getIntensPrecL(const);
%{
fit1.adjustWeightL_byPhase(matrix,10,false);
fit1.adjustWeightL_byRefL([[0 0 0]],0.1,true);
res=proj.fit();
%}
%res=proj.sweep();
%{
figure
contour(res.yL,res.xL,res.X2_redM')
set(gca,'xscale','log')
%}

%res=proj.vectorSweep();
%proj.setSweepVal(iSweepBestL);
%res=proj.fit();

%proj.sweepFit();
toc
%% Calculations


%{
peak1Info.H_refBragg=-[1;-1;3]/2;
peak1Info.xAng=0/180*pi;%additional rotation in X-Y plane to match experiment
%peak1Info.zoneAng=-0.13;%Tilt from zone axis
peak1Info.excAng=0;%0.Additional tilt from Bragg condition
peak1Info.iBragg=0;

Ug1=getUg(const,sim1,peak1Info.H_refBragg);
Xsi1=1/(diff1.lambda*Ug1);

%
peak2Info.H_refBragg=2*[6;-6;-4];
Ug2=getUg(const,sim1,peak2Info.H_refBragg);
Xsi2=1/(diff1.lambda*Ug2);

orientDoubleBragg(const,sim1,peak1Info,peak2Info);

g1BraggV=sim1.region.B_refM'*peak1Info.H_refBragg;
g2BraggV=sim1.region.B_refM'*peak2Info.H_refBragg;

sg1=getExcErr(g1BraggV,sim1.diff.kV,sim1.diff.nFoilV);
sg2=getExcErr(g2BraggV,sim1.diff.kV,sim1.diff.nFoilV);

[sg1 Xsi1]
[sg2 Xsi2]

sim1.H_refL=[];
sim1.H_refL=[0 0 0;peak1Info.H_refBragg';peak2Info.H_refBragg'];

%Thickness scan

scanInfo.nThick=200;
scanInfo.thickMin=0.1;
scanInfo.thickMax=500;
[thickResL,I_resL]=thicknessScan(const,basis,sim1,peak1Info,scanInfo);

%Rocking curve
scanInfo.nTilt=201;
scanInfo.tiltMin=-0.016;
scanInfo.tiltMax=0.016;
[tiltResL,I_resL]=rockingCurve(const,basis,sim1,peak1Info,scanInfo);
%}

%%Single SAED plot
%{
region1.orient(const,diff1);
region1.adjust();
sim1.expandBeamL(const);
sim1.prepare(const);
sim1.getIntensL(const);
nBeam=size(sim1.intensL,1);
sim1.createPlotList(1:nBeam,false,diff1.gMax);

diff1.autoExp=false;
sim1.plot(const);
%}
%%Single precession pattern

sim1.getIntensPrecL(const);
fit1.scaleIntens(const)
sim1.plot(const);

