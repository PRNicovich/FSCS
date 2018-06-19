function FSCSSimulation(particleBehavior, NRuns, crossTalk, BoxX, BoxY, tmax, fNameStart, varargin)

% 2D membrane diffusion simulations for FSCS manuscript
% 6 emission channels given by Lo and Ld reference spectra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs

% Ref spectra
% refLoSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\LUVsB_PSM_Chol_NR12S_2_spectrum.txt';
% refLdSpectrum = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\SLBs_DOPC_Chol_NR12S_3_spectrum.txt';

refLoSpectrum = '.\ReferenceSpectra\LUVsB_PSM_Chol_NR12S_2_spectrum.txt';
refLdSpectrum = '.\ReferenceSpectra\SLBs_DOPC_Chol_NR12S_3_spectrum.txt';


% Behavior of individual particles
% 'independent' - each particle has its own identity (Lo or Ld) and does not change
% 'staticDomains' - particle position in or out of domain determines
%                   identity, but domains do not move
% 'motileDomains' - particle position in or out of domain determines
%                   identity and domains translate randomly.  Particles
%                   themselves are static within a domain.
% 'switchingDomains' - particles switch between Lo and Ld while remaining
%                   static.  No spatial info in switching events.
% 'switchingDomainsWithIntermediate' - particles switch between Lo and Ld while remaining
%                   static.  No spatial info in switching events.
%                   Domains pass through 'intermediate' phase where equal
%                   Lo and Ld character.
% 'switchingGaussianDomains' - regions switch between Lo and Ld while dyes
%                   remain static. Switching done by appearance and
%                   disappearance of randomly-distributed Gaussian-shaped
%                   domains.
% 'staticBlendedMask' - static mask, with blurring parameter to blend
%                   binary mask into pseudo-continuous one 
% particleBehavior = 'switchingDomains';

    switch particleBehavior
        
        case 'independent'
            LoDiffCoeff = varargin{1};
            LdDiffCoeff = varargin{2};
        case 'staticDomains'
            LoDiffCoeff = varargin{1};
            LdDiffCoeff = varargin{2};
            IsingParams = varargin{3};
        case 'motileDomains'
            motileDomainDiffCoeff = varargin{1};
            IsingParams = varargin{2};
        case 'switchingDomains'
            LoToLd = varargin{1};
            LdToLo = varargin{2};
            diffFlag = varargin{3};
            if diffFlag == 1
                LoDiffCoeff = varargin{4};
                LdDiffCoeff = varargin{5};
            end
        case 'switchingDomainsWithIntermediate'
            LoToLd = varargin{1};
            LdToLo = varargin{2};
            IntermediateOutRate = varargin{3};
            diffFlag = varargin{4};
            if diffFlag == 1
                LoDiffCoeff = varargin{5};
                LdDiffCoeff = varargin{6};
            end
            
        case 'staticBlendedMask'
            LoDiffCoeff = varargin{1};
            LdDiffCoeff = varargin{2};
            IsingParams = varargin{3};
            GaussianSigma = varargin{4};
            
        case 'staticErodedMask'
            LoDiffCoeff = varargin{1};
            LdDiffCoeff = varargin{2};
            IsingParams = varargin{3};
            erosionMagnitude = varargin{4};
            
        case 'motileDomainsWithDiffusion'

            LoDiffCoeff = varargin{1};
            LdDiffCoeff = varargin{2};
            motileDomainDiffCoeff = varargin{3};
            IsingParams = varargin{4};
            
            
    end


% NRuns = 1; % Number of times to repeat simulation

% Include cross-talk between spectral channels if desired
% [LdInLo LoInLd]
% crossTalk = [.1 .1];

%Time units = ms
% tmax = 1e5;
tstep = 1; %  Timestep - best to keep this at 1
corrMax = 10e3; % Value to truncate correlation curve (ms)

% Simulation box size in um
% BoxX = 2;
% BoxY = 2;

%Gaussian detection volume
%Aspect controlls radius of spot between center and exp(-2) of gaussian beam shape.
%Aspect in um
aspectX = 0.230;
aspectY = 0.230;

Amp = 100; % Intensity amplitude (arb units)

%Particles per entire FOV
Nparts = 200;
       

% For 'staticDomains' or 'switchingDomains', a 2D Ising spin model is used
% Ising region is centered at (0,0) 
% IsingParams.InteractionStrength = 1.3806e-23; % Interaction strength in J
% IsingParams.Temperature = 0.5; % Kelvin
% IsingParams.PreRunIterations = 50; % tstep; Number of frames to evolve before starting FSCS simulation
% IsingParams.RegionSize = [999, 999]; % nm; Size of Ising box to simulate with random outside

%Diffusion coefficient in um^2/sec

% Diffusion steps to allow before experiment begins
% Used if particleBehavior == 'staticDomains'
DiffusionEvolutionSteps = 10000;
% From http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2966005/
% diffusion of lipids is ~1-3 um^2/sec
% LoDiffCoeff = 3;
% LdDiffCoeff = 1;

% Used if particleBehavior == 'motileDomains'
% from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2599850/
% Vesicle near membrane diff coefficient 5-8e-8 cm^2/sec
% motileDomainDiffCoeff = 0.5; 

% Used if particleBehavior == 'switchingDomains'
% Domain switching rates
% All in units of 1/sec
% rateConstant = ln(2)/t_{1/2}
% LoToLd = 5;
% LdToLo = 15;
% diffFlag = 1;

% Used if particleBehavior == 'switchingDomainsWithIntermediate'
% Above rates define transition rate out of Lo or Ld domain.  This
% describes rate from intermediate domain to either Lo or Ld (equal chance
% of each). 
% IntermediateOutRate = 0.1;

DiffSlowDownFactor = 0.01; % unused

% Used if particleBehavior == 'switchingGaussDomains'
% Lo domains in a background of Ld is simulated as a bunch of 
% 3D Gaussians in XYT space.  Each time step moves through T dimension.
% Lo defined as region with value above 1/(e^2).
% NLoDomains = 1e5;
% % Domain switching rates
% % Domain radius in um
% domRad = 0.05;
% % Domain length in ms (time units)
% domLength = 10;

% Plot colors
LoCurve = [46, 204, 90]/255;
LdCurve = [142, 68, 173]/255;
CrossCorr = [34, 126, 230]/255;

% LoMarker = 'o';
% LdMarker = 's';
% CrossMarker = 'd';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Lo, Ld ref spectra
LoData = dlmread(refLoSpectrum, '\t', 4, 0); 
Lo = LoData(:,2);
LdData = dlmread(refLdSpectrum, '\t', 4, 0);
Ld = LdData(:,2);

% Correct colors for plotting
% Reduce saturation and lighten all
LoLowSatCurve = hsv2rgb(rgb2hsv(LoCurve).*[1 0.4 1.2]);
LdLowSatCurve = hsv2rgb(rgb2hsv(LdCurve).*[1 0.4 1.3]);
CCLowSatCurve = hsv2rgb(rgb2hsv(CrossCorr).*[1 0.4 1.1]);

coorKeep = zeros(round(min([corrMax tmax])/tstep), 4, NRuns);

for N = 1:NRuns

    ParticleLocation = zeros(tmax, Nparts, 2);


    %Initial posistions are randomly assigned inside the box in X and Y
    ParticleLocation(1, :, 1) = rand(Nparts, 1)*BoxX - BoxX/2;
    ParticleLocation(1, :, 2) = rand(Nparts, 1)*BoxY - BoxY/2;

    lowboundX = -BoxX/2; 
    highboundX = BoxX/2;
    lowboundY = -BoxY/2; 
    highboundY = BoxY/2;

    BoxArea = (BoxX*BoxY)*1e-12; % in m^2


    particleIdentity = zeros(tmax, Nparts);

    % Set particle behavior 
    % Pure Lo is a value equal to 0, pure Ld equal to 1
    % Values between show fraction of Lo vs Ld character
    switch particleBehavior
        case 'independent'
            % Each particle assigned an identity which it keeps for entire
            % experiment
            particleIdentity = repmat(round(rand(1, Nparts)), tmax, 1);

        case 'staticDomains'

            IsingParams.Nimages = 1;
            IsingImgStack = IsingSimulationStack(IsingParams);
            particleMask = round(rand(BoxX*1000 + 1, BoxY*1000 + 1));
            particleMask((BoxX/2e-3) + (ceil(-IsingParams.RegionSize(1)/2):floor(IsingParams.RegionSize(1)/2)), ...
                (BoxY/2e-3) + (ceil(-IsingParams.RegionSize(2)/2):floor(IsingParams.RegionSize(2)/2))) = IsingImgStack;

            % Initialize particleIdentity
            particleIdentity = zeros(tmax, Nparts);
            for p = 1:Nparts
                particleIdentity(1, p) = particleMask(round(1000*(ParticleLocation(1, p, 1) + BoxX/2) + 1), ...
                    round(1000*(ParticleLocation(1, p, 2) + BoxY/2)) + 1);
            end


        case 'motileDomains'
            % Make a single big image that translates between time steps

            IsingParams.Nimages = 1;
            IsingImgStack = IsingSimulationStack(IsingParams);
            particleMask = round(rand(BoxX*1000 + 1, BoxY*1000 + 1));
            particleMask((BoxX/2e-3) + (ceil(-IsingParams.RegionSize(1)/2):floor(IsingParams.RegionSize(1)/2)), ...
                (BoxY/2e-3) + (ceil(-IsingParams.RegionSize(2)/2):floor(IsingParams.RegionSize(2)/2))) = IsingImgStack;

            % Initialize particleIdentity
            % Value doens't change through entire acquisition, so this can be
            % done for all frames at once.
            particleIdentity = zeros(tmax, Nparts);
            for p = 1:Nparts
                particleIdentity(:, p) = particleMask(round(1000*(ParticleLocation(1, p, 1) + BoxX/2) + 1), ...
                    round(1000*(ParticleLocation(1, p, 2) + BoxY/2)) + 1);
            end

        case 'switchingDomains'
            % Particles can switch between Lo and Ld by Hidden Markov Model

            RateConstants = [(1 - LoToLd/(tstep*1000)) LoToLd/(tstep*1000); 
                             LdToLo/(tstep*1000) (1 - LdToLo/(tstep*1000))];


            Emissions = [0 1;
                         1 0];

            Symbols = [0 1];

            particleIdentity = zeros(tmax, Nparts);

            for m = 1:Nparts

                [seq, states] = hmmgenerate(tmax/tstep, RateConstants, Emissions, 'symbols', Symbols);
                particleIdentity(:,m) = seq;

            end
            
        case 'switchingDomainsWithIntermediate'
            
            RateConstants = [(1 - LoToLd/(tstep*1000)) LoToLd/(tstep*1000) 0; 
                             0.5*IntermediateOutRate/(tstep*1000) (1 - IntermediateOutRate/(tstep*1000)) 0.5*IntermediateOutRate/(tstep*1000);
                            0 LdToLo/(tstep*1000) (1 - LdToLo/(tstep*1000))];


            Emissions = [0 0 1;
                         0 1 0;
                         1 0 0];

            Symbols = [0 0.5 1];

            particleIdentity = zeros(tmax, Nparts);

            for m = 1:Nparts

                [seq, states] = hmmgenerate((tmax + DiffusionEvolutionSteps)/tstep, RateConstants, Emissions, 'symbols', Symbols);
                particleIdentity(:,m) = 0.5*(states((DiffusionEvolutionSteps + 1):end)-1);

            end

        case 'switchingGaussianDomains'

            % Randomly choose XYT centers for NLoDomains domains
            xCent = BoxX*rand(NLoDomains, 1) - BoxX/2;
            yCent = BoxY*rand(NLoDomains, 1) - BoxY/2;
            tCent = (tmax/tstep)*rand(NLoDomains, 1) - (tmax/(tstep*2));

            % Particle identity is determined by presence in one of these
            % domains modeled as a XYT Gaussian

            particleIdentity(1,:) = double(any((1/(((2*pi)^(2/3))*(domRad*domRad*domLength))) * exp(-(((repmat(ParticleLocation(1,:,1), NLoDomains, 1) - repmat(xCent, 1, Nparts)).^2)/(2*domRad.^2) + ...
                ((repmat(ParticleLocation(1,:,2), NLoDomains, 1) - repmat(yCent, 1, Nparts)).^2)/(2*domRad^2) + ...
                (((1*ones(NLoDomains, Nparts)) - repmat(tCent, 1, Nparts)).^2)/(2*domLength^2))) > exp(-2)));
            
        case 'staticBlendedMask'

            IsingParams.Nimages = 1;
            IsingImgStack = IsingSimulationStack(IsingParams);

            particleMask = round(rand(BoxX*1000 + 1, BoxY*1000 + 1));
            particleMask((BoxX/2e-3) + (ceil(-IsingParams.RegionSize(1)/2):floor(IsingParams.RegionSize(1)/2)), ...
                (BoxY/2e-3) + (ceil(-IsingParams.RegionSize(2)/2):floor(IsingParams.RegionSize(2)/2))) = IsingImgStack;
            
            
            % Blur mask by GaussianSigma
            blurFilt = fspecial('gaussian', 2*GaussianSigma + 1, GaussianSigma);
            particleMask = imfilter(particleMask, blurFilt, 'symmetric', 'same');

            % Initialize particleIdentity
            particleIdentity = zeros(tmax, Nparts);
            for p = 1:Nparts
                particleIdentity(1, p) = particleMask(round(1000*(ParticleLocation(1, p, 1) + BoxX/2) + 1), ...
                    round(1000*(ParticleLocation(1, p, 2) + BoxY/2)) + 1);
            end
            
        case 'staticErodedMask'

            IsingParams.Nimages = 1;
            IsingImgStack = IsingSimulationStack(IsingParams);

            particleMask = round(rand(BoxX*1000 + 1, BoxY*1000 + 1));
            particleMask((BoxX/2e-3) + (ceil(-IsingParams.RegionSize(1)/2):floor(IsingParams.RegionSize(1)/2)), ...
                (BoxY/2e-3) + (ceil(-IsingParams.RegionSize(2)/2):floor(IsingParams.RegionSize(2)/2))) = IsingImgStack;
            
            
            % Erode mask by erosionMagnitude
            particleMask = bwmorph(particleMask, 'erode', erosionMagnitude);

            % Initialize particleIdentity
            particleIdentity = zeros(tmax, Nparts);
            for p = 1:Nparts
                particleIdentity(1, p) = particleMask(round(1000*(ParticleLocation(1, p, 1) + BoxX/2) + 1), ...
                    round(1000*(ParticleLocation(1, p, 2) + BoxY/2)) + 1);
            end
            
        case 'motileDomainsWithDiffusion'
            % Make a single big image that translates between time steps

            IsingParams.Nimages = 1;
            IsingImgStack = IsingSimulationStack(IsingParams);
            particleMask = round(rand(BoxX*1000 + 1, BoxY*1000 + 1));
            particleMask((BoxX/2e-3) + (ceil(-IsingParams.RegionSize(1)/2):floor(IsingParams.RegionSize(1)/2)), ...
                (BoxY/2e-3) + (ceil(-IsingParams.RegionSize(2)/2):floor(IsingParams.RegionSize(2)/2))) = IsingImgStack;

            % Distance of domain motion from starting point
            domainMotionCumulative = zeros(tmax, 2);
            prL = zeros(tmax, 2);
            
            % Initialize particleIdentity
            % Value doens't change through entire acquisition, so this can be
            % done for all frames at once.
            particleIdentity = zeros(tmax, Nparts);
            for p = 1:Nparts
                particleIdentity(:, p) = particleMask(round(1000*(ParticleLocation(1, p, 1) + BoxX/2) + 1), ...
                    round(1000*(ParticleLocation(1, p, 2) + BoxY/2)) + 1);
            end

    end

    % Add in cross talk between channels
    particleIdentity = particleIdentity + (1 - particleIdentity)*crossTalk(1) - particleIdentity*crossTalk(2);

    if strcmp(particleBehavior, 'staticDomains') | strcmp(particleBehavior, 'staticBlendedMask') | ...
            strcmp(particleBehavior, 'staticErodedMask') | strcmp(particleBehavior, 'staticErodedMask') | ...
            strcmp(particleBehavior, 'motileDomainsWithDiffusion');
        % Evolve diffusion parameters to capture any 'trapping' that happens
        % due to difference in diffusion time in different domains
        % repeatedly update particleIdentity and ParticleLocation in position 1
        % on each iteration of the loop
        for t = 1:DiffusionEvolutionSteps

            for p = 1:Nparts
                particleIdentity(1, p) = particleMask(round(1000*(ParticleLocation((1), p, 1) + BoxX/2) + 1), ...
                    round(1000*(ParticleLocation((1), p, 2) + BoxY/2)) + 1);
            end

            D = particleIdentity(1,:)*LoDiffCoeff + (1 - particleIdentity(1,:))*LdDiffCoeff;

            thetaDir = rand(Nparts, 1)*(2*pi) - pi;

            % rms displacement <r^2> is related to diffusion constant by <r^2> = 4 * D * \delta T
            % http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2966005/ eqn 5
            % um = um + sqrt(((um^2/sec) * (ms) * (1 sec/1000 ms)))
            ParticleLocation(1, :, 1) = ParticleLocation(1, :, 1) + sqrt((4*D*(tstep/1000))).*cos(thetaDir(:))';
            ParticleLocation(1, :, 2) = ParticleLocation(1, :, 2) + sqrt((4*D*(tstep/1000))).*sin(thetaDir(:))';

            ParticleLocation(1, (ParticleLocation(1, :, 1) > highboundX), 1) = ...
                ParticleLocation(1, (ParticleLocation(1, :, 1) > highboundX), 1) - BoxX;

            ParticleLocation(1, (ParticleLocation(1, :, 1) < lowboundX), 1) = ...
                ParticleLocation(1, (ParticleLocation(1, :, 1) < lowboundX), 1) + BoxX;

            ParticleLocation(1, (ParticleLocation(1, :, 2) > highboundY), 2) = ...
                ParticleLocation(1, (ParticleLocation(1, :, 2) > highboundY), 2) - BoxY;

            ParticleLocation(1, (ParticleLocation(1, :, 2) < lowboundY), 2) = ...
                ParticleLocation(1, (ParticleLocation(1, :, 2) < lowboundY), 2) + BoxY;

        end

    end
    

    %Calculating particle trajectories
    for t = 2:tstep:tmax



        % Update particleIdentity, if required
        switch particleBehavior

            case 'independent'
                % Do nothing.  Identity stays constant

            case 'staticDomains'

                % Particle identity in this frame is dependent on position in
                % the particleMask (to nearest nm) in previous frame
                for p = 1:Nparts
                    particleIdentity(t, p) = particleMask(round(1000*(ParticleLocation((t-1), p, 1) + BoxX/2) + 1), ...
                        round(1000*(ParticleLocation((t-1), p, 2) + BoxY/2)) + 1);
                end

            case 'motileDomains'
                % Identity stays constant

            case 'switchingDomains'
                % Identity already assigned

            case 'switchingGaussianDomains'

                particleIdentity(t,:) = double(any((1/(((2*pi)^(2/3))*(domRad*domRad*domLength))) * exp(-(((repmat(ParticleLocation(t,:,1), NLoDomains, 1) - repmat(xCent, 1, Nparts)).^2)/(2*domRad.^2) + ...
                    ((repmat(ParticleLocation(t,:,2), NLoDomains, 1) - repmat(yCent, 1, Nparts)).^2)/(2*domRad^2) + ...
                    (((t*ones(NLoDomains, Nparts)) - repmat(tCent, 1, Nparts)).^2)/(2*domLength^2))) > exp(-2)));
                
            case 'staticBlendedMask'

                % Particle identity in this frame is dependent on position in
                % the particleMask (to nearest nm) in previous frame
                for p = 1:Nparts
                    particleIdentity(t, p) = particleMask(round(1000*(ParticleLocation((t-1), p, 1) + BoxX/2) + 1), ...
                        round(1000*(ParticleLocation((t-1), p, 2) + BoxY/2)) + 1);
                end
                
            case 'staticErodedMask'

                % Particle identity in this frame is dependent on position in
                % the particleMask (to nearest nm) in previous frame
                for p = 1:Nparts
                    particleIdentity(t, p) = particleMask(round(1000*(ParticleLocation((t-1), p, 1) + BoxX/2) + 1), ...
                        round(1000*(ParticleLocation((t-1), p, 2) + BoxY/2)) + 1);
                end
                
            case 'motileDomainsWithDiffusion'
                % Identity changes with diffusion of particles
                try
                    for p = 1:Nparts
                        
                        % Update particle location based on position of
                        % particles with domain motion removed
                        partRefLocation = [round(1000*(ParticleLocation((t-1), p, 1) - domainMotionCumulative((t-1), 1) + BoxX/2) + 1), ...
                            round((1000*(ParticleLocation((t-1), p, 2) - domainMotionCumulative((t-1), 2) + BoxY/2)) + 1)];
                        
                        if partRefLocation(1) > 1000*BoxX
                            partRefLocation(1) = partRefLocation(1) - BoxX*1000;
                        end

                        if partRefLocation(1) < 1
                            partRefLocation(1) = partRefLocation(1) + BoxX*1000;
                        end

                        if partRefLocation(2) > 1000*BoxY
                            partRefLocation(2) = partRefLocation(2) - BoxY*1000;
                        end

                        if partRefLocation(2) < 1
                            partRefLocation(2) = partRefLocation(2) + BoxY*1000;
                        end
                        
                        prL(t, :) = partRefLocation;
                        
                        particleIdentity(t, p) = particleMask(partRefLocation(1), partRefLocation(2));
                        
                    end
                catch mError
                    assignin('base', 'p', p);
                    assignin('base', 't', t);
                    assignin('base', 'domainMotionCumulative', domainMotionCumulative);
                    assignin('base', 'particleMask', particleMask);
                    assignin('base', 'particleIdentity', particleIdentity);
                    assignin('base', 'ParticleLocation', ParticleLocation);

                    rethrow(mError);
                end

        end


        

    %     if strcmp(particleBehavior, 'switchingDomains')
    %         D = D*DiffSlowDownFactor;
    %     end

        if strcmp(particleBehavior, 'motileDomains')
            % In case of motile domains, the domain motion drives particle
            % motion

            thetaDir = rand(1, 1)*(2*pi) - pi;


            % Domain motion is in nm/timeStep
            domainMotion = round(1000*sqrt((4*motileDomainDiffCoeff*tstep/1000)).*[cos(thetaDir), sin(thetaDir)]);
            

            % Back to um
            ParticleLocation(t, :, 1) = ParticleLocation((t-1), :, 1) + round(domainMotion(1))/1000;
            ParticleLocation(t, :, 2) = ParticleLocation((t-1), :, 2) + round(domainMotion(2))/1000;

    
        elseif strcmp(particleBehavior, 'motileDomainsWithDiffusion')
            
            thetaDir = rand(1, 1)*(2*pi) - pi;

            % Domain motion is in nm/timeStep
            domainMotion = round(1000*sqrt((4*motileDomainDiffCoeff*tstep/1000)).*[cos(thetaDir), sin(thetaDir)]);
            
            domainMotionCumulative(t, :) = domainMotionCumulative((t-1), :) + round(domainMotion)/1000; % In um, range [-Box/2 Box/2]
            
            if domainMotionCumulative(t, 1) > BoxX/2
                domainMotionCumulative(t, 1) = domainMotionCumulative(t, 1) - BoxX;
            end
            
            if domainMotionCumulative(t, 1) < -BoxX/2
                domainMotionCumulative(t, 1) = domainMotionCumulative(t, 1) + BoxX;
            end
            
            if domainMotionCumulative(t, 2) > BoxY/2
                domainMotionCumulative(t, 2) = domainMotionCumulative(t, 2) - BoxY;
            end
            
            if domainMotionCumulative(t, 2) < -BoxY/2
                domainMotionCumulative(t, 2) = domainMotionCumulative(t, 2) + BoxY;
            end
            
            
            
            
            % Particle motion due to diffusion
            D = particleIdentity(t,:)*LoDiffCoeff + (1 - particleIdentity(t,:))*LdDiffCoeff;
            
            thetaDir = rand(Nparts, 1)*(2*pi) - pi;

            % rms displacement <r^2> is related to diffusion constant by <r^2> = 4 * D * \delta T
            % http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2966005/ eqn 5
            % um = um + sqrt(((um^2/sec) * (ms) * (1 sec/1000 ms)))
            ParticleLocation(t, :, 1) = ParticleLocation((t-1), :, 1) + sqrt((4*D*(tstep/1000))).*cos(thetaDir(:))' + round(domainMotion(1))/1000;
            ParticleLocation(t, :, 2) = ParticleLocation((t-1), :, 2) + sqrt((4*D*(tstep/1000))).*sin(thetaDir(:))' + round(domainMotion(2))/1000;


        elseif strcmp(particleBehavior, 'switchingGaussianDomains')

            % Don't move as domains switch around particle
            
        elseif strcmp(particleBehavior, 'switchingDomains') & diffFlag == 0;

            % Don't move as domains switch around particle
            
        elseif strcmp(particleBehavior, 'switchingDomainsWithIntermediate')  & diffFlag == 0;

            % Don't move as domains switch around particle

        else

            % staticDomains or independent or staticBlendedMask or
            % staticErodedMask

            D = particleIdentity(t,:)*LoDiffCoeff + (1 - particleIdentity(t,:))*LdDiffCoeff;
            
            thetaDir = rand(Nparts, 1)*(2*pi) - pi;

            % rms displacement <r^2> is related to diffusion constant by <r^2> = 4 * D * \delta T
            % http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2966005/ eqn 5
            % um = um + sqrt(((um^2/sec) * (ms) * (1 sec/1000 ms)))
            ParticleLocation(t, :, 1) = ParticleLocation((t-1), :, 1) + sqrt((4*D*(tstep/1000))).*cos(thetaDir(:))';
            ParticleLocation(t, :, 2) = ParticleLocation((t-1), :, 2) + sqrt((4*D*(tstep/1000))).*sin(thetaDir(:))';

        end

        ParticleLocation(t, (ParticleLocation(t, :, 1) > highboundX), 1) = ...
            ParticleLocation(t, (ParticleLocation(t, :, 1) > highboundX), 1) - BoxX;

        ParticleLocation(t, (ParticleLocation(t, :, 1) < lowboundX), 1) = ...
            ParticleLocation(t, (ParticleLocation(t, :, 1) < lowboundX), 1) + BoxX;

        ParticleLocation(t, (ParticleLocation(t, :, 2) > highboundY), 2) = ...
            ParticleLocation(t, (ParticleLocation(t, :, 2) > highboundY), 2) - BoxY;

        ParticleLocation(t, (ParticleLocation(t, :, 2) < lowboundY), 2) = ...
            ParticleLocation(t, (ParticleLocation(t, :, 2) < lowboundY), 2) + BoxY;


        if mod((10*t), tmax) == 0;
            NComp = (100*t/tmax);
            fprintf(1,'Calculating particles: %d %%',NComp);
            fprintf('\n');
        end

    end


    % Determine flourescence intensity at each time point for total number of
    % particles

    intByFluor = Amp * exp(-(((ParticleLocation(:,:,1).^2)/(2*aspectX^2)) + ...
        ((ParticleLocation(:,:,2).^2)/(2*aspectY^2))));


    LoFluor = intByFluor.*particleIdentity;
    LdFluor = intByFluor.*(1 - particleIdentity);


    % Convert this to TTTR for processing
    spectralFluor = sum(LoFluor, 2)*Lo' + sum(LdFluor, 2)*Ld';

    corrOut = FSCSintensityCorrelation(spectralFluor, Lo, Ld);

    corrOut(corrOut(:,1) > corrMax, :) = [];
    if le(corrMax, tmax)
        corrOut(end, :, :) = [];
    end

    assignin('base', 'coorKeep', coorKeep);
    assignin('base', 'corrOut', corrOut);
    
    coorKeep(:,:,N) = corrOut;
    
   
    switch particleBehavior
        
        case 'independent'

        case 'staticDomains'
            
            % Save mask image
            fileName = sprintf('%s%s_%s_N=%d.tif', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, N);
            imwrite(particleMask, fileName, 'tiff');

        case 'motileDomains'

        case 'switchingDomains'

        case 'switchingGaussianDomains'
            
        case 'staticBlendedMask'
            
            % Save mask image
            fileName = sprintf('%s%s_%s_N=%d.tif', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, N);
            imwrite(particleMask, fileName, 'tiff');
            
        case 'staticErodedMask'
            
            % Save mask image
            fileName = sprintf('%s%s_%s_N=%d.tif', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, N);
            imwrite(particleMask, fileName, 'tiff');
            
        case 'motileDomainsWithDiffusion'
            
            % Save mask image
            fileName = sprintf('%s%s_%s_N=%d.tif', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, N);
            imwrite(particleMask, fileName, 'tiff');
        
    end

end

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
switch particleBehavior
    case 'independent'
        
    case 'staticDomains'
        
    case 'motileDomains'
        
    case 'switchingDomains'
        
    case 'switchingGaussianDomains'
        
    case 'staticBlendedMask'
        
	case 'staticErodedMask'
        
    case 'motileDomainsWithDiffusion'
        
        assignin('base', 'domainMotionCumulative', domainMotionCumulative);
        assignin('base', 'prL', prL);
        
end

            

plotFig = figure();
plotAx = axes('parent', plotFig);
semilogx(plotAx, coorKeep(:, 1, 1), squeeze((coorKeep(:,2,:))), 'color', LoLowSatCurve, 'linewidth', 1);
hold on
semilogx(plotAx, coorKeep(:, 1, 1), squeeze((coorKeep(:,3,:))), 'color', LdLowSatCurve, 'linewidth', 1);
semilogx(plotAx, coorKeep(:, 1, 1), squeeze((coorKeep(:,4,:))), 'color', CCLowSatCurve, 'linewidth', 1);
LoMeanPlot = semilogx(plotAx, coorKeep(:, 1, 1), squeeze(mean(coorKeep(:,2,:), 3)), 'color', LoCurve, 'linewidth', 2);
LdMeanPlot = semilogx(plotAx, coorKeep(:, 1, 1), squeeze(mean(coorKeep(:,3,:), 3)), 'color', LdCurve, 'linewidth', 2);
CCMeanPlot = semilogx(plotAx, coorKeep(:, 1, 1), squeeze(mean(coorKeep(:,4,:), 3)), 'color', CrossCorr, 'linewidth', 2);
hold off
legend([LoMeanPlot, LdMeanPlot, CCMeanPlot], 'Lo', 'Ld', 'Cross');
set(plotAx, 'fontsize', 12);
set(plotFig, 'color', [1 1 1]);
xlabel(plotAx, '\tau (ms)', 'fontsize', 12);
ylabel(plotAx, 'Correlation', 'fontsize', 12);

if sum(crossTalk(:)) > 0
    particleBehavior = sprintf('%s_ctalk_%.2f,%.2f', particleBehavior, crossTalk(1), crossTalk(2));
end


% Save figure
fileName = sprintf('%s%s_%s_N%d.tif', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, NRuns);
print(plotFig, '-dtiff', fileName);
fileName = sprintf('%s%s_%s_N%d.svg', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, NRuns);
print(plotFig, '-dsvg', fileName);


close(plotFig);

% Save curves w/ metadata
fileName = sprintf('%s%s_%s_N%d.txt', 'C:\Users\Rusty\Documents\MATLAB\FSCS\FSCS\', fNameStart, particleBehavior, NRuns);
fID = fopen(fileName, 'w+');
fprintf(fID, '# ParticleBehavior : %s\r\n', particleBehavior);
fprintf(fID, '# NRuns : %d\r\n', NRuns);
fprintf(fID, '# Crosstalk : %.2f %.2f\r\n', crossTalk(1), crossTalk(2));
fprintf(fID, '# Tmax : %d\r\n', tmax);
fprintf(fID, '# Tstep : %d\r\n', tstep);
fprintf(fID, '# corrMax : %.0f\r\n', corrMax);
fprintf(fID, '# BoxX BoxY : %.1f %.1f\r\n', BoxX, BoxY);
fprintf(fID, '# FocalVolumeSize : %.4f x %.4f\r\n', aspectX, aspectY);
fprintf(fID, '# Intensity : %.0f\r\n', Amp);
fprintf(fID, '# Nparticles : %df\r\n', Nparts);
fprintf(fID, '# DiffusionEvolutionSteps : %d\r\n', DiffusionEvolutionSteps);

switch particleBehavior
    case 'independent'
        fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff);
    case 'staticDomains'
        fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff); 
        fprintf(fID, '# Ising IntStrength Temp PreRun RegionSizeX RegionSizeY : %.5f %.2f %d %d %d\r\n', ...
            IsingParams.InteractionStrength, IsingParams.Temperature, IsingParams.PreRunIterations, IsingParams.RegionSize);
    case 'motileDomains'
        fprintf(fID, '# MotileDomainCoeff : %.2f\r\n', motileDomainDiffCoeff);
        fprintf(fID, '# Ising IntStrength Temp PreRun RegionSizeX RegionSizeY : %.5f %.2f %d %d %d\r\n', ...
            IsingParams.InteractionStrength, IsingParams.Temperature, IsingParams.PreRunIterations, IsingParams.RegionSize);
    case 'switchingDomains'
        fprintf(fID, '# TransitionCoeff LoLd LdLo : %.2f %.2f\r\n', LoToLd, LdToLo);
        fprintf(fID, '# DiffFlag : %.d\r\n', diffFlag); 
        if diffFlag == 1
            fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff);  
        end
    case 'switchingDomainsWithIntermediate'
        fprintf(fID, '# TransitionCoeff LoLd Intermed LdLo : %.2f %.2f %.2f\r\n', LoToLd, LdToLo);
        fprintf(fID, '# DiffFlag : %.d\r\n', diffFlag); 
        if diffFlag == 1
            fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff);  
        end
        
     case 'staticBlendedMask'
        fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff); 
        fprintf(fID, '# Ising IntStrength Temp PreRun RegionSizeX RegionSizeY : %.5f %.2f %d %d %d\r\n', ...
            IsingParams.InteractionStrength, IsingParams.Temperature, IsingParams.PreRunIterations, IsingParams.RegionSize);
        fprintf(fID, '# GaussianSigma : %.2f\r\n', GaussianSigma);
        
	case 'staticErodedMask'
        fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff); 
        fprintf(fID, '# Ising IntStrength Temp PreRun RegionSizeX RegionSizeY : %.5f %.2f %d %d %d\r\n', ...
            IsingParams.InteractionStrength, IsingParams.Temperature, IsingParams.PreRunIterations, IsingParams.RegionSize);
        fprintf(fID, '# Erosion Magnitude : %.2f\r\n', erosionMagnitude);
        
	case 'motileDomainsWithDiffusion'
        fprintf(fID, '# DiffCoeff Lo Ld : %.2f %.2f\r\n', LoDiffCoeff, LdDiffCoeff); 
        fprintf(fID, '# Ising IntStrength Temp PreRun RegionSizeX RegionSizeY : %.5f %.2f %d %d %d\r\n', ...
            IsingParams.InteractionStrength, IsingParams.Temperature, IsingParams.PreRunIterations, IsingParams.RegionSize);
        fprintf(fID, '# MotileDomainCoeff : %.2f\r\n', motileDomainDiffCoeff);
end
% fprintf(fID, '# NLoDomains domRad domLength : %d %.4f %.2f\r\n', NLoDomains, domRad, domLength);

for k = 1:size(coorKeep, 3)
    
    fprintf(fID, '####################################\r\n');
    fprintf(fID, 'Time (ms)\tLoAuto\tLdAuto\tCrossCorr\r\n');
    
    for m = 1:size(coorKeep, 1)
        
        fprintf(fID, '%.1f\t%.6f\t%.6f\t%.6f\r\n', squeeze(coorKeep(m, :, k)));
        
    end
end

fclose(fID);

end

% save Fluor.txt Fluor -ascii;
% clear all;

%plot(time,AllPartsY{1},'b.',time,AllPartsX{1},'g.');
%axis([0 tmax -BoxY/2 BoxY/2]);



