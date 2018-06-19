
% Behavior of individual particles
% 'independent' - each particle has its own identity (Lo or Ld) and does not change
% 'staticDomains' - particle position in or out of domain determines
%                   identity, but domains do not move
% 'motileDomains' - particle position in or out of domain determines
%                   identity and domains translate randomly.  Particles
%                   themselves are static within a domain.
% 'switchingDomains' - particles switch between Lo and Ld while remaining
%                   static.  No spatial info in switching events.
% 'switchingGaussianDomains' - regions switch between Lo and Ld while dyes
%                   remain static. Switching done by appearance and
%                   disappearance of randomly-distributed Gaussian-shaped
%                   domains.

particleBehavior = 'independent';

NRuns = 20; % Number of times to repeat simulation

% Include cross-talk between spectral channels if desired
% [LdInLo LoInLd]
crossTalk = [0 0];

tmax = 1e5;

IsingParams.InteractionStrength = 1.3806e-23; % Interaction strength in J
IsingParams.Temperature = 0.5; % Kelvin
IsingParams.PreRunIterations = 50; % tstep; Number of frames to evolve before starting FSCS simulation
IsingParams.RegionSize = [1999, 1999]; % nm; Size of Ising box to simulate with random outside

BoxX = 5;
BoxY = 5;

LoDiffCoeff = 0.1;
LdDiffCoeff = 0.3;

LoToLd = 5;
LdToLo = 50;

% Used if particleBehavior == 'motileDomains'
% from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2599850/
% Vesicle near membrane diff coefficient 5-8e-8 cm^2/sec
motileDomainDiffCoeff = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start calling conditions
% 
% % Independent diffusion, no switching, plus crossTalk
fprintf(1, 'FSCS sim independent no cross-talk\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_NoCtalk', LoDiffCoeff, LdDiffCoeff);

fprintf(1, 'FSCS sim independent low cross-talk\n');
crossTalk = [0.05 0.05];
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_LowCTalk', LoDiffCoeff, LdDiffCoeff);

fprintf(1, 'FSCS sim independent medium cross-talk\n');
crossTalk = [0.1 0.1];
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MidCTalk', LoDiffCoeff, LdDiffCoeff);

fprintf(1, 'FSCS sim independent high cross-talk\n');
crossTalk = [0.2 0.2];
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_HighCTalk', LoDiffCoeff, LdDiffCoeff);
% % 
% % Static domains w/ different diffusion parameters
particleBehavior = 'staticDomains'; 
crossTalk = [0 0];
LoDiffCoeff = 0.3;
LdDiffCoeff = 0.3;
fprintf(1, 'FSCS sim static equal diffusion rates\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_EqualDiff', LoDiffCoeff, LdDiffCoeff, IsingParams);

LoDiffCoeff = 0.1;
LdDiffCoeff = 0.3;
fprintf(1, 'FSCS sim static Ld faster diffusion\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_LdFaster', LoDiffCoeff, LdDiffCoeff, IsingParams);

LoDiffCoeff = 0.3;
LdDiffCoeff = 0.1;
fprintf(1, 'FSCS sim static Lo faster diffusion\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_LoFaster', LoDiffCoeff, LdDiffCoeff, IsingParams);

% Static domains, different domain sizes
Be sure to double-check the resulting domain sizes with these
LoDiffCoeff = 0.1;
LdDiffCoeff = 0.3;
IsingParams.Temperature = 0.001; % Kelvin
IsingParams.PreRunIterations = 100;
fprintf(1, 'FSCS sim static Extreme domains\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_ExtremeDomains', LoDiffCoeff, LdDiffCoeff, IsingParams);

IsingParams.Temperature = 0.01; % Kelvin
fprintf(1, 'FSCS sim static Big domains\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_BigDomains', LoDiffCoeff, LdDiffCoeff, IsingParams);

IsingParams.Temperature = 50; % Kelvin
IsingParams.PreRunIterations = 50;
fprintf(1, 'FSCS sim static Small domains\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SmallDomains', LoDiffCoeff, LdDiffCoeff, IsingParams);

% Static domains, throw in some crosstalk
IsingParams.Temperature = 0.5; % Kelvin
crossTalk = [0.1 0.1];
fprintf(1, 'FSCS sim static with crosstalk\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_CrossTalkDomains', LoDiffCoeff, LdDiffCoeff, IsingParams);

% 
% % 14082016
% % Motile domains starting point
particleBehavior = 'motileDomains';
crossTalk = [0 0];
LoDiffCoeff = 0.01;
LdDiffCoeff = 0.01; 
motileDomainDiffCoeff = 1; 
IsingParams.Temperature = 0.5; % Kelvin
IsingParams.PreRunIterations = 50;
IsingParams.Nimages = 1;
IsingParams.RegionSize = [999, 999]; % nm; Size of Ising box to simulate with random outside


% % Medium domains, different diffusion speeds
fprintf(1, 'FSCS sim motile domains medium diffusion\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MedDomMidDiff', motileDomainDiffCoeff, IsingParams);

motileDomainDiffCoeff = 5;
fprintf(1, 'FSCS sim motile domains fast diffusion\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MedDomFastDiff', motileDomainDiffCoeff, IsingParams);

motileDomainDiffCoeff = 1/5;
fprintf(1, 'FSCS sim  motile domains slow diffusion\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MedDomSlowDiff', motileDomainDiffCoeff, IsingParams);

% Small domains
IsingParams.Temperature = 50; % Kelvin
motileDomainDiffCoeff = 1; 
fprintf(1, 'FSCS sim motile domains small domains\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SmallDomMidDiff', motileDomainDiffCoeff, IsingParams);

% Big domains
IsingParams.Temperature = 0.01; % Kelvin
IsingParams.PreRunIterations = 100;
fprintf(1, 'FSCS sim motile domains big domains\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_BigDomMidDiff', motileDomainDiffCoeff, IsingParams);

% Extreme domains
IsingParams.Temperature = 0.001; % Kelvin
fprintf(1, 'FSCS sim motile domains extreme domains\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_ExtremeDomMidDiff', motileDomainDiffCoeff, IsingParams);

% 
%%%%%%
% Switching domains
particleBehavior = 'switchingDomains';
LoToLd = 5;
LdToLo = 15;
diffFlag = 1;
IntermediateOutRate = 1;

fprintf(1, 'FSCS sim switching domains Lo slower\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchLoSlower', LoToLd, LdToLo, diffFlag, LoDiffCoeff, LdDiffCoeff);

LoToLd = 15;
LdToLo = 15;
fprintf(1, 'FSCS sim switching domains equal rates\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchEqualRates', LoToLd, LdToLo, diffFlag, LoDiffCoeff, LdDiffCoeff);

LoToLd = 15;
LdToLo = 5;
fprintf(1, 'FSCS sim switching domains Ld slower\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchLdSlower', LoToLd, LdToLo, diffFlag, LoDiffCoeff, LdDiffCoeff);

LoToLd = 5;
LdToLo = 15;
diffFlag = 0;
fprintf(1, 'FSCS sim switching domains Lo slower with no diffusion\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchLoSlowerNoDiffusion', LoToLd, LdToLo, diffFlag, LoDiffCoeff, LdDiffCoeff);

%%%%%%
% Switching domains with intermediates
particleBehavior = 'switchingDomainsWithIntermediate';
LoToLd = 5;
LdToLo = 15;
IntermediateOutRate = 1;
diffFlag = 1;

fprintf(1, 'FSCS sim switching domains with intermediate, slow out\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchIntermedSlow', LoToLd, LdToLo, IntermediateOutRate, diffFlag, LoDiffCoeff, LdDiffCoeff);

IntermediateOutRate = 10;
fprintf(1, 'FSCS sim switching domains with intermediate, mid out\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchIntermedMid', LoToLd, LdToLo, IntermediateOutRate, diffFlag, LoDiffCoeff, LdDiffCoeff);

IntermediateOutRate = 45;
fprintf(1, 'FSCS sim switching domains with intermediate, fast out\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_SwitchIntermedFast', LoToLd, LdToLo, IntermediateOutRate, diffFlag, LoDiffCoeff, LdDiffCoeff);


%%%%%%%%%%%%%%%%%%%%55
% More sims after looking at previous batch of results
% Switching domains and static domains look about right with peak in
% cross-correlation, but the amplitude is incorrect.  Need to find a way to
% get cross-corr lower than both Lo and Ld.  

% Static domains
% Lots of cross-talk between channels
particleBehavior = 'staticDomains';
crossTalk = [0.2 0.2];
fprintf(1, 'FSCS sim static with 20% crosstalk\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MedCrossTalkDomains', LoDiffCoeff, LdDiffCoeff, IsingParams);

crossTalk = [0.4 0.4];
fprintf(1, 'FSCS sim static with 20% crosstalk\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_HighCrossTalkDomains', LoDiffCoeff, LdDiffCoeff, IsingParams);

% Maybe binary phase mask is incorrect for this approach.  Introducing new
% phaseBehavior, 'staticBlendedMask', with parameters for blurring of Ising
% mask with Gaussian kernel.
crossTalk = [0 0];
particleBehavior = 'staticBlendedMask';
GaussianSigma = 1;
fprintf(1, 'FSCS sim static blended\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Blended1px', LoDiffCoeff, LdDiffCoeff, IsingParams, GaussianSigma);

GaussianSigma = 5;
fprintf(1, 'FSCS sim static blended\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Blended5px', LoDiffCoeff, LdDiffCoeff, IsingParams, GaussianSigma);

GaussianSigma = 10;
fprintf(1, 'FSCS sim static blended\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Blended10px', LoDiffCoeff, LdDiffCoeff, IsingParams, GaussianSigma);

GaussianSigma = 10;
fprintf(1, 'FSCS sim static blended\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Blended10pxa', LoDiffCoeff, LdDiffCoeff, IsingParams, GaussianSigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All of these have had equal proportion of two phases within image area.
% Trying an erosion step to see if this changes anything 
particleBehavior = 'staticErodedMask';
erosionMagnitude = 1;
fprintf(1, 'FSCS sim static eroded by 1 pixel\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Eroded1px', LoDiffCoeff, LdDiffCoeff, IsingParams, erosionMagnitude);

particleBehavior = 'staticErodedMask';
erosionMagnitude = 5;
fprintf(1, 'FSCS sim static eroded by 5 pixels\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Eroded5px', LoDiffCoeff, LdDiffCoeff, IsingParams, erosionMagnitude);

particleBehavior = 'staticErodedMask';
erosionMagnitude = 10;
fprintf(1, 'FSCS sim static eroded by 10 pixels\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_Eroded10px', LoDiffCoeff, LdDiffCoeff, IsingParams, erosionMagnitude);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combining motile domains with dye diffusion to see if this gets us the
% right relationship between auto-correlation and cross-correlation
% amplitudes
particleBehavior = 'motileDomainsWithDiffusion';

NRuns = 20;

crossTalk = [0 0];
LoDiffCoeff = 0.01;
LdDiffCoeff = 0.03; 
motileDomainDiffCoeff = 0.5; 
IsingParams.Temperature = 0.001; % Kelvin % Big domain temp and pre-runs
IsingParams.PreRunIterations = 100;
IsingParams.Nimages = 1;
IsingParams.RegionSize = [1999, 1999]; % nm; Size of Ising box to simulate with random outside

% % Medium domains, different diffusion speeds
fprintf(1, 'FSCS sim motile domains and diffusion particles\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MotileAndDiffusionHandChosenValues', LoDiffCoeff, LdDiffCoeff, motileDomainDiffCoeff, IsingParams);

% Turns out this works.  With the right combo of motileDomainDiffCoeff and
% LoDiffCoeff, LdDiffCoeff you can get the right order of amplitudes
LoDiffCoeff = 0.1;
LdDiffCoeff = 0.3; 
motileDomainDiffCoeff = 1; 
fprintf(1, 'FSCS sim motile domains and diffusion particles standard\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MotileAndDiffusionStandDiff', LoDiffCoeff, LdDiffCoeff, motileDomainDiffCoeff, IsingParams);

LoDiffCoeff = 0.01;
LdDiffCoeff = 0.03; 
motileDomainDiffCoeff = 0.1; 
fprintf(1, 'FSCS sim motile domains and diffusion particles 10x slower\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MotileAndDiffusion10xSlower', LoDiffCoeff, LdDiffCoeff, motileDomainDiffCoeff, IsingParams);

LoDiffCoeff = 0.001;
LdDiffCoeff = 0.003; 
motileDomainDiffCoeff = 0.1; 
fprintf(1, 'FSCS sim motile domains and diffusion particles slower diff\n');
FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MotileAndDiffusionSlowerDiff', LoDiffCoeff, LdDiffCoeff, motileDomainDiffCoeff, IsingParams);


% particleBehavior = 'motileDomains';
% fprintf(1, 'FSCS sim motile domains medium diffusion\n');
% FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, 'FSCSSim_MedDomMidDiff', motileDomainDiffCoeff, IsingParams);
