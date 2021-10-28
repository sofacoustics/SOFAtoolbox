% miro :: Measured Impulse Response Object
%
% miro Class Definition
%
% MIRO V1, Release 1.04 - 07/Feb/2014
% MIRO V1, Release 1.05 - 25/Dez/2016
% MIRO V1, Release 1.06 - 20/Sep/2017

% #Author: Benjamin Bernschütz (2012-2014)
% #Author: Johannes M. Arend (2014-2017)
% #Author: Christoph Pörschmann (TH Köln, University of Applied Sciences,
%   Institute of Communications Engineering, 
%   Mail christoph.poerschmann@th-koeln.de, Phone +49 221 8275 2495) (2012-2017)

% #Author: Chris Pike, BBC R&D (02.2014)
% #Author: Johannes M. Arend, TH Köln (12.2016)
% #Author: Johannes M. Arend & Christoph Pörschmann, TH Köln (09.2017)
% #Author: Michael Mihocic: header documentation updated (28.10.2021)

classdef miro
    
    
    properties
        
        name               = [];
        context            = [];
        location           = [];
        date               = [];
        engineer           = [];
        contact            = [];
        comments           = [];
        miroVersion        = [];
        type               = [];
        fs                 = [];
        taps               = [];
        nIr                = [];
        excitationSignal   = [];
        gapTime            = [];
        microphone         = [];
        source             = [];
        audioInterface     = [];
        micPreamp          = [];
        capturingSystem    = [];
        systemLoopLatency  = [];
        latencyCompensated = [];
        headCut            = [];
        sourcePosition     = [];
        sourceDistance     = [];
        avgAirTemp         = [];
        avgRelHumidity     = [];
        positionReference  = [];
        postProcessing     = [];
        normalization      = [];
        ctIRnormalization  = [];
        quadGrid           = [];
        scatterer          = [];
        radius             = [];
        azimuth            = [];
        elevation          = [];
        quadWeight         = [];
        chOne              = [];
        chTwo              = [];
        irChOne            = [];
        irChTwo            = [];
        centerMicrophone   = [];
        irCenter           = [];
        returnTaps         = [];
        resampleToFS       = [];
        headWin            = [];
        tailWin            = [];
        headPhone          = [];
        hpcfKernel         = [];
        headPhoneComp      = false;
        mic                = [];
        mcfKernel          = [];
        micComp            = false;
        shutUp             = false;
        angles             = 'RAD';
        
    end
    
    methods
        
        function [ir, azimuth, elevation, quadWeight] = getIR(obj, irID)
            
            % [ir, azimuth, elevation, quadWeight] = getIR(obj, irID)
            %
            % Returns the impulse response(s), angles and the
            % quadrature weight for a specific ID number [irID].
            % The returned impulse response has two channels
            % if the object contains two channels.
            %
            % irID = 0 adresses the omni center IR. The center IR is
            %          not filtered by the headphone filter.
            
            if ~isempty(obj.resampleToFS)
                resampleIR = true;
                gComDiv = gcd(obj.fs, obj.resampleToFS);
                p = obj.resampleToFS / gComDiv;
                q = obj.fs / gComDiv;
            else
                resampleIR = false;
            end
            
            if irID == 0
                
                irOne      = obj.irCenter(:,1);
                irOne      = cast(irOne, 'double');
                azimuth    = [];
                elevation  = [];
                quadWeight = [];
                
            else
                
                irOne      = obj.irChOne(:,irID);
                irOne      = cast(irOne, 'double');
                azimuth    = obj.azimuth(irID);
                elevation  = obj.elevation(irID);
                quadWeight = obj.quadWeight(irID);
            end
            
            if obj.headPhoneComp && irID > 0      % Applying HP Filter
                if size(obj.hpcfKernel,2)==1
                    irOne = miro_fastConv(irOne, obj.hpcfKernel);
                else
                    irOne = miro_fastConv(irOne, obj.hpcfKernel(:,1));
                end
            end
            
            if obj.micComp && irID > 0      % Applying Mic Filter %JMA-Edit
                irOne = miro_fastConv(irOne, obj.mcfKernel);
            end
            
            if obj.headWin > 0                    % Windowing Head
                irOne = miro_winHead(irOne, obj.headWin);
            end
            
            if obj.returnTaps < obj.taps          % Truncation
                irOne = irOne(1:obj.returnTaps);
            end
            
            if obj.tailWin > 0                    % Windowing Tail
                irOne = miro_winTail(irOne, obj.tailWin);
            end
            
            ir(:,1) = irOne;
            
            if ~isempty(obj.irChTwo) && irID > 0  % Channel Two
                
                irTwo = obj.irChTwo(:,irID);
                irTwo = cast(irTwo, 'double');
                
                if obj.headPhoneComp
                    if size(obj.hpcfKernel,2)==1
                        irTwo = miro_fastConv(irTwo, obj.hpcfKernel);
                    else
                        irTwo = miro_fastConv(irTwo, obj.hpcfKernel(:,2));
                    end
                end
                
                if obj.micComp                  % Applying Mic Filter %JMA-Edit
                    irTwo = miro_fastConv(irTwo, obj.mcfKernel);
                end
                
                if obj.headWin > 0                % Windowing Head
                    irTwo = miro_winHead(irTwo, obj.headWin);
                end
                
                if obj.returnTaps < obj.taps      % Truncation
                    irTwo = irTwo(1:obj.returnTaps);
                end
                
                if obj.tailWin > 0                % Windowing Tail
                    irTwo = miro_winTail(irTwo, obj.tailWin);
                end
                
                ir(:,2) = irTwo;
            end
            
            if resampleIR                         % Resampling to new FS
                ir = resample(ir,double(p),double(q));
            end
            
        end
        
        function [irID, azimuth, elevation] = closestIr(obj, az_approx, el_approx)
            
            % [irID, azimuth, elevation] = closestIr(obj, az_approx, el_approx)
            %
            % This method returns the ID number [irID] of the closest
            % angle to [az_approx, el_approx] that is available
            % within the object. [azimuth] and [elevation] return
            % the corresponding closest fitting angles
            
            
            if strcmp(obj.angles,'DEG')
                az_approx   = az_approx*pi/180;
                el_approx   = el_approx*pi/180;
                obj.shutUp  = true;
                objRAD      = setRAD(obj);
                quadrature  = getQuadrature(objRAD);
            else
                quadrature  = getQuadrature(obj);
            end
            
            [Xg,Yg,Zg]   = sph2cart(quadrature(:,1),quadrature(:,2)-pi/2,1);
            [Xp,Yp,Zp]   = sph2cart(az_approx,el_approx-pi/2,1);
            distances    = sqrt((Xg-Xp).^2+(Yg-Yp).^2+(Zg-Zp).^2);
            irID         = find(distances == min(distances));
            azimuth      = obj.azimuth(irID);
            elevation    = obj.elevation(irID);
            
        end
        
        function quadrature = getQuadrature(obj)
            
            % quadrature = getQuadrature(obj)
            %
            % Returns the quadrature including azimuth and elevation
            % angles and weights.
            
            quadrature(:,1) = double(obj.azimuth);
            quadrature(:,2) = double(obj.elevation);
            quadrature(:,3) = double(obj.quadWeight);
            
        end
        
        function plotQuadrature(obj)
            
            % plotQuad(obj)
            %
            % Shows a sphere plot of the quadrature
            
            obj.shutUp = true;
            obj        = setRAD(obj);
            quadrature = getQuadrature(obj);
            
            [x,y,z] = sph2cart(quadrature(:,1),quadrature(:,2)-pi/2,1.01);
            
            colormap Gray;
            
            if size(x,1)>1500
                plot3(x,y,z,'marker','.','markerfacecolor','g','color','g','linestyle','none')
            else
                plot3(x,y,z,'marker','o','markerfacecolor','g','color','g','linestyle','none')
            end
            
            axis off; hold on; grid off; sphere; axis equal; rotate3d on; light;
            alpha(.6); lighting phong; hold off; set(gcf,'color','w');
            
            if strcmp(obj.type,'HRIR')
                title(['Sampling Points or Virtural Source Positions in: ',strrep(obj.name,'_','-')]);
            elseif strcmp(obj.type,'BRIR')
                title(['Head orientations in: ',strrep(obj.name,'_','-')]); %JMA-Edit: Distinguish between HRIR and BRIR type
            end
        end
        
        function miroCoordinates(obj)
            
            % miroCoordinates(obj)
            %
            % This method illustrates the coordinate system
            % that is used by miro
            
            RAD = true;
            
            if nargin > 0
                if strcmp(obj.angles,'DEG')
                    RAD = false;
                end
            end
            
            delete(gca)
            
            line([-5,5],[0,0],[0,0],'Color','r')
            hold on
            line([0,0],[-5,5],[0,0],'Color','g')
            line([0,0],[0,0],[-5,5],'Color','b')
            
            if RAD
                front  = 'FRONT: AZ 0/2pi - EL pi/2';
                back   = 'BACK: AZ pi - EL pi/2';
                left   = 'LEFT: AZ pi/2 - EL pi/2';
                right  = 'RIGHT: AZ 3pi/2 - EL pi/2';
                top    = 'TOP: AZ x - EL 0';
                bottom = 'BOTTOM: AZ x - EL pi/2';
            else
                front  = 'FRONT: AZ 0/360° - EL 90°';
                back   = 'BACK: AZ 180° - EL 90°';
                left   = 'LEFT: AZ 90° - EL 90°';
                right  = 'RIGHT: AZ 270° - EL 90°';
                top    = 'TOP: AZ x - EL 0°';
                bottom = 'BOTTOM: AZ x - EL 180°';
            end
            
            text([7,7],  [0,0],  [0,0],  front, 'HorizontalAlignment','center')
            text([-7,-7],[0,0],  [0,0],  back,  'HorizontalAlignment','center')
            text([0,0],  [7,7],  [0,0],  left,  'HorizontalAlignment','center')
            text([0,0],  [-7,-7],[0,0],  right, 'HorizontalAlignment','center')
            text([0,0],  [0,0],  [7,7],  top,   'HorizontalAlignment','center')
            text([0,0],  [0,0],  [-7,-7],bottom,'HorizontalAlignment','center')
            
            set(gca,'Projection','perspective'); set(gcf,'color','w');
            box off; axis off; view(-60,34);
            title('miro ::: Coordinate System');
        end
        
        function obj = setDEG(obj)
            
            % obj = setDEG(obj)
            %
            % Changes the object's angle reference to DEG.
            % All angle handling is then in DEG.
            %
            % ! The default angle reference is RAD.
            
            if strcmp(obj.angles,'RAD');
                obj.angles    = 'DEG';
                obj.azimuth   = obj.azimuth*180/pi;
                obj.elevation = obj.elevation*180/pi;
                if ~obj.shutUp
                    fprintf('Angle reference is now DEG.\n');
                end
            else
                if ~obj.shutUp
                    fprintf('Angle reference is DEG. (No Change)\n');
                end
            end
        end
        
        function obj = setRAD(obj) % Set Angles to RAD
            
            % obj = setRAD(obj)
            %
            % Changes the object's angle reference to RAD.
            % All angle handling is then in RAD.
            %
            % ! The default angle reference is RAD.
            
            if strcmp(obj.angles,'DEG');
                obj.angles    = 'RAD';
                obj.azimuth   = obj.azimuth*pi/180;
                obj.elevation = obj.elevation*pi/180;
                if ~obj.shutUp
                    fprintf('Angle reference is now RAD.\n');
                end
            else
                if ~obj.shutUp
                    fprintf('Angle reference is RAD. (No Change)\n');
                end
            end
        end
        
        function obj = setReturnTaps(obj, returnTaps, tailWin)
            
            % obj = setReturnTaps(obj, returnTaps, tailWin)
            %
            % Allows for changing the number of returned impulse
            % response taps. By default, all available taps are
            % returned (retunedTaps = taps). Setting the returnedTaps
            % property to a smaller value will cause that all returned
            % impulse responses are cut off at the value defined by
            % [returnTaps] and are by default windowed using a
            % half-sided Hann window of the size (returnTaps/8).
            % The size of the window can be defined by [tailWin].
            % Using [tailWin = 0] turns off windowing.
            
            
            if returnTaps <= obj.taps
                obj.returnTaps = returnTaps;
                
                if nargin < 3
                    obj.tailWin = round(returnTaps/8);
                else
                    obj.tailWin = tailWin;
                    if tailWin> obj.returnTaps
                        obj.tailWin = obj.returnTaps;
                        if ~obj.shutUp
                            fprintf('Warning: The window cannot be larger than returnTaps!\n')
                        end
                    end
                end
                if ~obj.shutUp
                    fprintf('returnTaps and tailWin set\n');
                end
            else
                if ~obj.shutUp
                    fprintf('ERROR: returnTaps value too large.\n');
                end
            end
        end
        
        function obj = setResampling(obj, targetFS)
            
            % obj = setResampling(obj, targetFS)
            %
            % Sets the resampleToFS property and can be used to
            % extract the IRs at a different audio sampling
            % rate. The resampling has effect on the methods:
            % getIR(), playAudio() and dropWaveFile().
            %
            % targetFS defines the target sampling rate, e.g. 44100Hz
            %
            % If called without the targteFS argument, the resampleToFS
            % property is set to [] (default) and disabled.
            %
            % The resampling process is done within the getIR method
            % after truncation and windowing. The returned IRs then do
            % NOT have the amount of samples set in returnTaps property
            % as these refer to the original FS.
            %
            % The resampling is based on the native MATLAB resample()
            % method which is included in the Signal Processing Library.
            %
            % ADVICE: If not urgently necessary try to avoid resampling.
            
            if nargin < 2
                
                obj.resampleToFS = [];
                if ~obj.shutUp
                    fprintf(['Resampling disabled.\n']);
                end
                
            else
                obj.resampleToFS = round(targetFS);
                if ~obj.shutUp
                    fprintf(['Resampling enabled and new target FS is set.\n']);
                end
            end
        end
        
        function dropWaveFile(obj, irID, nBits, filename)
            
            % dropWaveFile(obj, irID, nBits, filename)
            %
            % Drops a .WAV audio file containing the impulse
            % responses for a specific ID number [irID].
            % The arguments bitdepht [nBits] and filename
            % [filename] are optional.
            %
            % Defaults:
            % ------------------------------------------------
            % nBits    = 16
            % filename = obj.name,'_IR',num2str(irID), ...
            %            '-AZ',az, 'EL',num2str(el), ...
            %            obj.angles, [obj.headphone]
            % ------------------------------------------------
            % WARNING: NO DITHERING/NOISESHAPING APPLIED
            
            if nargin < 3
                nBits = 16;
            end
            
            if nargin < 4
                if strcmp(obj.angles,'RAD')
                    az = num2str(round(obj.azimuth(irID)*1e4));
                    el = num2str(round(obj.elevation(irID)*1e4));
                else
                    az = num2str(round(obj.azimuth(irID)*1e2));
                    el = num2str(round(obj.elevation(irID)*1e2));
                end
                filename = [obj.name,'_IR',num2str(irID),'-AZ',az,'EL',num2str(el),obj.angles];
                
                if obj.headPhoneComp
                    filename = [filename,'_',obj.headPhone];
                end
                
                if obj.micComp
                    filename = [filename,'_',obj.mic]; %JMA Edit - MCF (needed for NarDasS research project)
                end
            end
            
            if strcmp(filename,'NoCoor') %JMA Edit - Small hack to export wav files without any info on coordinates in the filename (needed for NarDasS research project)
                filename = ['SSR_',obj.name];
                
                if obj.headPhoneComp
                    filename = [filename,'_',obj.headPhone];
                end
                
                if obj.micComp
                    filename = [filename,'_',obj.mic]; %JMA Edit - MCF (needed for NarDasS research project)
                end
            end
                
            ir = getIR(obj, irID);
            
            if isempty(obj.resampleToFS)
                %JMA-Edit - Use wavwrite for matlab < 2012b (8.0) and audiowrite for matlab > 2012b
                if verLessThan ('matlab','8.0')
                    wavwrite(ir, obj.fs, nBits, filename); 
                else
                    audiowrite([filename,'.wav'], ir, obj.fs, 'BitsPerSample', nBits);
                end
            else
                %JMA-Edit - Use wavwrite for matlab < 2012b (8.0) and audiowrite for matlab > 2012b
                if verLessThan ('matlab','8.0')
                    wavwrite(ir, obj.resampleToFS, nBits, filename); 
                else
                    audiowrite([filename,'.wav'], ir, obj.resampleToFS, 'BitsPerSample', nBits);
                end
            end
            if ~obj.shutUp
                fprintf(['Wave file has been written.\n']);
            end
            
        end
        
        function playAudio(obj, irID, audioSignal)
            
            % playAudio(obj, irID, audioSignal)
            %
            % Plays an audio signal convolved with the
            % impulse response(s) given by irID.
            % The method serves for a quick pre-listening
            % of the datasets. If the signal [audioSignal]
            % has more than one channel, the first channel
            % is taken.
            
            if size(audioSignal,2) > size(audioSignal,1)
                audioSignal = audioSignal';
            end
            
            audioSignal = audioSignal(:,1);
            ir          = getIR(obj, irID);
            
            finalAudio(:,1) = miro_fastConv(audioSignal,ir(:,1));
            
            if size(ir,2) == 2
                finalAudio(:,2) = miro_fastConv(audioSignal,ir(:,2));
            end
            
            if ~obj.shutUp
                fprintf('Playing audio ...\n');
            end
            
            if isempty(obj.resampleToFS)
                soundsc(finalAudio, obj.fs);
            else
                soundsc(finalAudio, obj.resampleToFS);
            end
            
            if ~obj.shutUp
                fprintf('... done.\n');
            end
            
        end
        
        function obj = setHeadPhones(obj, hpcf, linearPhase)
            
            % obj = setHeadPhones(obj, hpcf, linearPhase)
            %
            % This method sets headphone compensation filters
            % for HRIR and BRIR datasets. The filters are
            % applied to all outgoing impulse responses while
            % obj.headphoneComp is set "true".
            
            if strcmp(obj.type,'BRIR') || strcmp(obj.type,'HRIR')
                
                if nargin < 3
                    linearPhase = false;
                end
                
                obj.headPhone = hpcf.hpName;
                
                if linearPhase
                    obj.hpcfKernel = hpcf.linPhase;
                else
                    obj.hpcfKernel = hpcf.minPhase; %DEFAULT
                end
                
                obj.headPhoneComp = true;
                if ~obj.shutUp
                    fprintf(['Headphone filters (',hpcf.hpName,') are set.\n']);
                end
            else
                if ~obj.shutUp
                    fprintf('Headphone filters can only be applied to objects of type HRIR or BRIR.');
                end
                return
                
            end
        end
        
        %JMA Edit - Microphone Compensation Filters (needed for NarDasS research project)
        function obj = setMic(obj, mcf, linearPhase) 
            
            % obj = setMic(obj, mcf, linearPhase)
            %
            % This method sets microphone compensation filters
            % for HRIR and BRIR datasets. The filters are
            % applied to all outgoing impulse responses while
            % obj.micComp is set "true".
            
            if strcmp(obj.type,'BRIR') || strcmp(obj.type,'HRIR')
                
                if nargin < 3
                    linearPhase = false;
                end
                
                obj.mic = mcf.micName;
                
                if linearPhase
                    obj.mcfKernel = mcf.linPhase;
                else
                    obj.mcfKernel = mcf.minPhase; %DEFAULT
                end
                
                obj.micComp = true;
                if ~obj.shutUp
                    fprintf(['Mic compensation filters (',mcf.micName,') are set.\n']);
                end
            else
                if ~obj.shutUp
                    fprintf('Mic compensation filters can only be applied to objects of type HRIR or BRIR.');
                end
                return
                
            end
        end
        
        function miroToSSR(obj, mirror, nBits, filename, normalize)
            
            % miroToSSR(obj, mirror, nBits, filename, normalize)
            %
            % This method writes a 720 channel interleaved
            % wave file for the Sound Scape Renderer (SSR):
            %
            % http://www.tu-berlin.de/?id=ssr
            %
            % The object must be a circular HRIR or BRIR set.
            %
            % Defaults:
            % ------------------------------------------------
            % mirror   = false (not mirrored)
            %            (This option can be useful for
            %            symmetrical venues to mirror a source
            %            e.g. from Left to Right.)
            % nBits    = 16
            % filename = 'SSR_', obj.name, [obj.headphone]
            % normalize= true
            % ------------------------------------------------
            % WARNING: NO DITHERING/NOISESHAPING APPLIED
            
            if (strcmp(obj.type,'BRIR')  || strcmp(obj.type,'HRIR')) && obj.nIr == 360
                
                if nargin < 2
                    mirror = false;
                end
                
                if nargin < 3
                    nBits = 16;
                end
                
                if nargin < 4
                    if mirror
                        filename = ['SSR_',obj.name,'_mirror'];
                    else
                        filename = ['SSR_',obj.name];
                    end
                    
                    if obj.headPhoneComp
                        filename = [filename,'_',obj.headPhone];
                    end
                    
                    if obj.micComp
                        filename = [filename,'_',obj.mic]; %JMA Edit -  Add mcf
                    end
                    
                end
                
                if nargin < 5
                    normalize = true;
                end
                
                if mirror
                    fprintf('Source mirrored ');
                end
                
                irArray = zeros(size(obj.getIR(1),1), 720);
                inc = 1;
                
                for i = 1:360
                    if ~mod(i,10)
                        fprintf('|');
                    end
                    
                    if mirror
                        IR = obj.getIR(361-i);
                        irArray(:,inc) = IR(:,2);
                        inc = inc+1;
                        irArray(:,inc) =  IR(:,1);
                        inc = inc+1;
                    else
                        IR = obj.getIR(i);
                        irArray(:,inc) = IR(:,1);
                        inc = inc+1;
                        irArray(:,inc) =  IR(:,2);
                        inc = inc+1;
                    end
                end
                
                fprintf('\n\n');
                
                if (normalize)
                    %disp('Normalizing...'); %JMA-Edit
                    irArray = 0.99*irArray/max(max(abs(irArray)));
                    disp('Normalization set true - peak normalization applied'); %JMA-Edit - Always print info about normalization
                else
                    disp('Normalization set false - no peak normalization applied');
                end
                
                miro_wavwrite([filename,'.wav'], irArray, obj.fs, nBits); %JMA-Edit - Use miro_wavwrite
                disp(['SSR file sucessfully generated: ', filename])
            else
                if ~obj.shutUp
                    fprintf('Failed: SSR export requires circular HRIR or BRIR objects.\n');
                end
                return
            end
        end
        
        function [timeDataCH1, timeDataCH2] = miroToSOFiA(obj)
            
            % timeData = miroToSOFiA(obj)
            %
            % Returns structs that are readable by the F/D/T
            % function of the SOFiA sound field analysis toolbox.
            % The following S/T/C transform core transforms
            % the object's data into the spherical harmonics domain.
            %
            % SOFiA: http://code.google.com/p/sofia-toolbox/
            
            obj = obj.setRAD;
            
            if ~isempty(obj.resampleToFS)
                downSampleFactor = obj.fs / obj.resampleToFS;
                fs = obj.resampleToFS;
            else
                downSampleFactor = 1;
                fs = obj.fs;
            end
            
            timeDataCH1.FS             = double(fs);
            timeDataCH1.radius         = double(obj.radius);
            timeDataCH1.quadratureGrid = obj.getQuadrature;
            timeDataCH1.downSample     = double(downSampleFactor);
            timeDataCH1.averageAirTemp = double(obj.avgAirTemp);
            timeDataCH1.irOverlay      = 0;
            timeDataCH1.centerIR       = obj.getIR(0)';
            
            timeDataCH1.impulseResponses = zeros(obj.nIr,size(timeDataCH1.centerIR,2));
            
            if ~isempty(obj.chTwo)
                twoChannels = true;
                timeDataCH2.FS             = fs;
                timeDataCH2.radius         = double(obj.radius);
                timeDataCH2.quadratureGrid = obj.getQuadrature;
                timeDataCH2.downSample     = double(downSampleFactor);
                timeDataCH2.averageAirTemp = double(obj.avgAirTemp);
                timeDataCH2.irOverlay      = 0;
                timeDataCH2.centerIR       = obj.getIR(0)';
                timeDataCH2.impulseResponses = zeros(obj.nIr,size(timeDataCH2.centerIR,2));
            else
                twoChannels = false;
                timeDataCH2 = [];
            end
            
            for irCnt = 1:obj.nIr
                
                impulseResponse = obj.getIR(irCnt);
                
                timeDataCH1.impulseResponses(irCnt,:) =  impulseResponse(:,1)';
                timeDataCH1.irOverlay = timeDataCH1.irOverlay + impulseResponse(:,1)';
                
                if twoChannels
                    timeDataCH2.impulseResponses(irCnt,:) =  impulseResponse(:,2)';
                end
            end
            
            timeDataCH1.irOverlay=sum(timeDataCH1.impulseResponses,1);
            timeDataCH1.irOverlay=abs(timeDataCH1.irOverlay/max(abs(timeDataCH1.irOverlay)));
            
            if twoChannels
                timeDataCH2.irOverlay=sum(timeDataCH2.impulseResponses,1);
                timeDataCH2.irOverlay=abs(timeDataCH2.irOverlay/max(abs(timeDataCH2.irOverlay)));
            end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Private methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ab = miro_fastConv(a,b)

% Internal use only

NFFT = size(a,1)+size(b,1)-1;
A    = fft(a,NFFT);
B    = fft(b,NFFT);
AB   = A.*B;
ab   = ifft(AB);

end

function ir = miro_winHead(ir, wLen)

% Internal use only

c = 1:wLen;
w = 0.5+0.5*cos(2*pi*(c-((2*wLen-1)/2))/(2*wLen-1));
ir(1:size(w,2)) = ir(1:size(w,2)) .* w';

end

function ir = miro_winTail(ir, wLen)

% Internal use only

c = wLen:2*wLen-1;
w = 0.5+0.5*cos(2*pi*(c-((2*wLen-1)/2))/(2*wLen-1));
ir(end-size(w,2)+1:end) = ir(end-size(w,2)+1:end) .* w';

end

%CP-Edit: wavwrite function to write multichannel wav files with more than
%512 channels. Required to write wavfiles for the SSR
function [] = miro_wavwrite(filename, Y, Fs, nbits)
%This function writes a multichannel .wav file
%Due to the removement of the wavread function it is currently not possible
%to write or read multichannel wav-files with extremely large number of
%channels (n > 512)
%To create BRIRs or HRIRs for the SSR (SoundScapeRenderer) multi-channel
%.wav-files are required
% Input Parameters: 
%       filename:   Filename of the .wav-File
%       Y:          Audio data of the .wav-file
%       Fs:         Sampling rate
%       nbits:      Number of bits per sample
% C.Pörschmann / TH Köln / January 2017


if nargin < 4    
    nbits=16;
end

% calculate som relevant values
number_of_channels=length(Y(1,:));
number_of_samples=length(Y(:,1));
FileSize=44 + nbits/8 * number_of_channels * number_of_samples;

%Open File
fileID = fopen(filename,'w');

% Write header
fwrite(fileID,'RIFF','*char');
fwrite(fileID,FileSize-8,'int');
fwrite(fileID,'WAVE','*char');
fwrite(fileID,'fmt ','*char');

FmtSize=16;
fwrite(fileID,FmtSize,'int');

Format=1;
fwrite(fileID,Format,'short');

fwrite(fileID,number_of_channels,'short');

fwrite(fileID,Fs,'int');

BytesPerSecond = Fs*nbits/8*number_of_channels;
fwrite(fileID,BytesPerSecond,'int'); 

Block_Alignment = number_of_channels * nbits/8;
fwrite(fileID,Block_Alignment,'short');


BitsPerSample = nbits ; 
fwrite(fileID,BitsPerSample,'short');

fwrite(fileID,'data','*char');

Datasize = FileSize - 44; 
fwrite(fileID,Datasize,'int');

% Change Amplitudes from [-1;1] as used in Matlab to integer
Y = Y * (2^(BitsPerSample-1)); 


for my_loop=1:size(Y,1)
    
    switch BitsPerSample
        case 8
            fwrite(fileID,Y(my_loop,:),'bit8');
        case 16
            fwrite(fileID,Y(my_loop,:),'bit16');
        case 24
            fwrite(fileID,Y(my_loop,:),'bit24');
        case 32
            fwrite(fileID,Y(my_loop,:),'bit32');
        otherwise
            error('Error writing SSR_file - Wrong number of bits per sample')
    end    
end
%Close File
fclose(fileID);

end

%CP-Edit: counterpart to miroWavread. Not really needed actually :)
function [Y, Fs, nbits] = miro_wavread(filename)
%This function reads a multichannel .wav file
%Due to the removement of the wavread function it is currently not possible
%to write or read multichannel wav-files with extremely large number of
%channels. 
%To create BRIRs or HRIRs for the SSR (SoundScapeRenderer) multi-channel
%.wav-files are required
% Input Parameters: 
%       filename:   Filename of the .wav-File
% Return Values: 
%       Y:          Audio data of the .wav-file
%       Fs:         Sampling rate
%       nbits:      Number of bits per sample

% C.Pörschmann / TH Köln / January 2017


%open File
fileID = fopen(filename,'r');

% read header
my_wave.header.RIFF= fread(fileID,4,'*char');
my_wave.header.RIFF=my_wave.header.RIFF';
my_wave.header.RiffSize = fread(fileID,1,'int');
my_wave.header.WAVE= fread(fileID,4,'*char');
my_wave.header.WAVE=my_wave.header.WAVE';
my_wave.header.FMT= fread(fileID,4,'*char');
my_wave.header.FMT=my_wave.header.FMT';
my_wave.header.FmtSize = fread(fileID,1,'int');
my_wave.header.Format = fread(fileID,1,'short');
my_wave.header.Channels = fread(fileID,1,'short');
my_wave.header.SamplesPerSecond = fread(fileID,1,'int');
my_wave.header.BytesPerSecond = fread(fileID,1,'int');
my_wave.header.BlockAlign = fread(fileID,1,'short');
my_wave.header.BitsPerSample = fread(fileID,1,'short');

%Adjust things, if header is longer than expected
my_wave.header.chunk = fread(fileID,my_wave.header.FmtSize-16,'*char');

my_wave.header.nextblock = fread(fileID,4,'*char');
my_wave.header.nextblock=my_wave.header.nextblock';

% IF a 'fact' chunk is included in the wavfile
if my_wave.header.nextblock=='fact'
    my_wave.header.fact_length = fread(fileID,1,'int');
    my_wave.header.chunk = fread(fileID,my_wave.header.fact_length,'*char');
    my_wave.header.nextblock = fread(fileID,4,'*char');
    my_wave.header.nextblock=my_wave.header.nextblock';
end
if my_wave.header.nextblock~='data'
    error ('This type of wav-format is not implemented');
end

my_wave.header.Datasize = fread(fileID,1,'int');
    
if my_wave.header.Format~=1
    error ('This type of wav-format is not implemented');
end


number_of_samples=my_wave.header.Datasize*8/my_wave.header.BitsPerSample/my_wave.header.Channels;
%read audio data
switch my_wave.header.BitsPerSample    
    case 8
        Y=fread(fileID,[my_wave.header.Channels number_of_samples],'bit8');
    case 16
        Y=fread(fileID,[my_wave.header.Channels number_of_samples],'bit16');
    case 24
        Y=fread(fileID,[my_wave.header.Channels number_of_samples],'bit24');
    case 32
        Y=fread(fileID,[my_wave.header.Channels number_of_samples],'bit32');
    otherwise
        error('Not implemented number of bits per sample in this wav-format')
end
Y=Y';

% Change amplitudes from integer to [-1;1] as used in Matlab 
Y = Y / (2^(my_wave.header.BitsPerSample-1)); 

%Calculate return values
Fs = my_wave.header.SamplesPerSecond;
nbits = my_wave.header.BitsPerSample;
%Close audio file
fclose(fileID);

end
