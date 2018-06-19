function [Channel, AbsoluteTime, MacroTime, MicroTime] =  pt3Import(filepath, varargin)

% MATLAB port of pt3import python method from :
% https://github.com/dwaithe/FCS_point_correlator/blob/master/import_methods.py
% Parallelization implemented to speed up processing.
% Can call this explicitly with 'Parallel' optional second argument (on by
% default) or use Picoquant implementation with 'PQ' as second argument.

% To do - implement editing out of overflow events from macrotime,
% microtime

if nargin == 2
    processFlag = varargin{1};
else
    processFlag = 'Parallel';
end

    f = fopen(filepath, 'rb');
    
    fseek(f, 0, 'eof');
    fileSize = ftell(f);
    fseek(f, 0, 'bof');
    
    
    Ident = fread(f, 16, 'uint8=>char')';
    FormatVersion = fread(f, 6, 'uint8=>char')';
    CreatorName = fread(f, 18, 'uint8=>char')';
    CreatorVersion = fread(f, 12, 'uint8=>char')';
    FileTime = fread(f, 18, 'uint8=>char')';
    CRLF = fread(f, 2, 'uint8=>char')';
    CommentField = fread(f, 256, 'uint8=>char')';
    
%     disp(ftell(f)) % 328 bytes
    
    Curves = fread(f, 1, 'int32');
    BitsPerRecord = fread(f, 1, 'int32');

    RoutingChannels = fread(f, 1, 'int32');
  
    NumberOfBoards = fread(f, 1, 'int32');

    ActiveCurve = fread(f, 1, 'int32');
 
    MeasurementMode = fread(f, 1, 'int32');

    SubMode = fread(f, 1, 'int32');

    RangeNo = fread(f, 1, 'int32');

    Offset = fread(f, 1, 'int32');

    AcquisitionTime = fread(f, 1, 'int32');

    StopAt = fread(f, 1, 'int32');
   
    StopOnOvfl = fread(f, 1, 'int32');
     
    Restart = fread(f, 1, 'int32');
          
    DispLinLog = fread(f, 1, 'int32');
         
    DispTimeFrom = fread(f, 1, 'int32');
        
    DispTimeTo = fread(f, 1, 'int32');
         
    DispCountFrom = fread(f, 1, 'int32');
           
    DispCountTo = fread(f, 1, 'int32');
            

%     disp(ftell(f)) % 72 more bytes (400)
             
    DispCurveMapTo = [];
    DispCurveShow =[];
    
    for i = 1:8 
        DispCurveMapTo = [DispCurveMapTo; fread(f, 1, 'int32')];
        
        DispCurveShow = [DispCurveShow; fread(f, 1, 'int32')];

    end
    
%     disp(ftell(f)) % 64 more bytes (464)
    
    ParamStart =[];
    ParamStep =[];
    ParamEnd =[];
    for i = 1:3; 
        ParamStart = [ParamStart; fread(f, 1, 'int32')];

        ParamStep = [ParamStep; fread(f, 1, 'int32')];
        
        ParamEnd = [ParamEnd; fread(f, 1, 'int32')];
       
    end
        
%     disp(ftell(f)) % 36 more bytes (500)
    
    RepeatMode = fread(f, 1, 'int32');
    RepeatsPerCurve = fread(f, 1, 'int32');
    RepeatTime = fread(f, 1, 'int32');
    RepeatWait = fread(f, 1, 'int32');
    ScriptName = fread(f, 20, 'uint8=>char')';

%     disp(ftell(f)) % 36 more bytes (536)
    
    % Board-specific header

    HardwareIdent = fread(f, 16, 'uint8=>char')';
    HardwareVersion = fread(f, 8, 'uint8=>char')';
    HardwareSerial = fread(f, 1, 'int32');
    SyncDivider = fread(f, 1, 'int32');

    CFDZeroCross0 = fread(f, 1, 'int32');
    CFDLevel0 = fread(f, 1, 'int32');
    CFDZeroCross1 = fread(f, 1, 'int32');
    CFDLevel1 = fread(f, 1, 'int32');

    Resolution = fread(f, 1, 'float32');
    
%     disp(ftell(f)) % 52 more bytes (588)

    % For format version 2.0

    RouterModelCode      = fread(f, 1, 'int32');
    RouterEnabled        = fread(f, 1, 'int32');

    % Router channel 1
    RtChan1_InputType    = fread(f, 1, 'int32');
    RtChan1_InputLevel   = fread(f, 1, 'int32');
    RtChan1_InputEdge    = fread(f, 1, 'int32');
    RtChan1_CFDPresent   = fread(f, 1, 'int32');
    RtChan1_CFDLevel     = fread(f, 1, 'int32');
    RtChan1_CFDZeroCross = fread(f, 1, 'int32');
    % Router channel 2
    RtChan2_InputType    = fread(f, 1, 'int32');
    RtChan2_InputLevel   = fread(f, 1, 'int32');
    RtChan2_InputEdge    = fread(f, 1, 'int32');
    RtChan2_CFDPresent   = fread(f, 1, 'int32');
    RtChan2_CFDLevel     = fread(f, 1, 'int32');
    RtChan2_CFDZeroCross = fread(f, 1, 'int32');
    % Router channel 3
    RtChan3_InputType    = fread(f, 1, 'int32');
    RtChan3_InputLevel   = fread(f, 1, 'int32');
    RtChan3_InputEdge    = fread(f, 1, 'int32');
    RtChan3_CFDPresent   = fread(f, 1, 'int32');
    RtChan3_CFDLevel     = fread(f, 1, 'int32');
    RtChan3_CFDZeroCross = fread(f, 1, 'int32');
    % Router channel 4
    RtChan4_InputType    = fread(f, 1, 'int32');
    RtChan4_InputLevel   = fread(f, 1, 'int32');
    RtChan4_InputEdge    = fread(f, 1, 'int32');
    RtChan4_CFDPresent   = fread(f, 1, 'int32');
    RtChan4_CFDLevel     = fread(f, 1, 'int32');
    RtChan4_CFDZeroCross = fread(f, 1, 'int32');

   % disp(ftell(f)) % 136 more bytes (724)
    
    % T3 mode specific header.

    ExtDevices = fread(f, 1, 'int32');

    Reserved1 = fread(f, 1, 'int32');
    Reserved2 = fread(f, 1, 'int32');
    CntRate0 = fread(f, 1, 'int32');
    CntRate1 = fread(f, 1, 'int32');

    StopAfter = fread(f, 1, 'int32');
    StopReason = fread(f, 1, 'int32');
    Records = fread(f, 1, 'uint32');
    ImgHdrSize = fread(f, 1, 'int32');
    
   % disp(ftell(f)) % 36 more bytes (760)

    % Special Header for imaging.
    if ImgHdrSize > 0
        ImgHdr = fread(f, ImgHdrSize, 'uint8=>char')';
    end
    ofltime = 0;

    cnt_1=0; cnt_2=0; cnt_3=0; cnt_4=0; cnt_Ofl=0; cnt_M=0; cnt_Err=0; % counters
    WRAPAROUND=65536;

    % Put file Save info here.

    syncperiod = 1e9/CntRate0;
    % outfile stuff here.
    % fpout.
    % T3RecordArr = [];
    
%     chanArr = zeros(Records, 1);
%     trueTimeArr = zeros(Records, 1);
%     dTimeArr= zeros(Records, 1);
    
%     disp(fileSize)
%     disp(ftell(f))
%     disp(Records)
    % f1=open('./testfile', 'w+')
%     for b = 1:Records
%         T3Record = fread(f, 1, 'uint32');
%         
%         % T3RecordArr.append(T3Record)
%         
%         nsync = bitand(uint32(T3Record), 65535);
%         
%         chan = bitand(uint32(bitsra(T3Record, 28)), 15);
%         
%         chanArr(b) = chan;
%         
%         % f1.write(str(i)+" "+str(T3Record)+" "+str(nsync)+" "+str(chan)+" ")
%         dtime = 0;
%         
%         if chan == 1
%             cnt_1 = cnt_1+1;
%             dtime = bitand(uint32(bitsra(T3Record, 16)), 4095);
%             % f1.write(str(dtime)+" ")
%         elseif chan == 2 
%             cnt_2 = cnt_2+1;
%             dtime = bitand(uint32(bitsra(T3Record, 16)), 4095);
%             %f1.write(str(dtime)+" ")
%         elseif chan == 3 
%             cnt_3 = cnt_3+1;
%             dtime = bitand(uint32(bitsra(T3Record, 16)), 4095);
%             %f1.write(str(dtime)+" ")
%         elseif chan == 4 
%             cnt_4 = cnt_4+1;
%             dtime = bitand(uint32(bitsra(T3Record, 16)), 4095);
%             %f1.write(str(dtime)+" ")
%         elseif chan == 15
%             markers = bitand(uint32(bitsra(T3Record, 16)), 4095);
%             
%             if markers == 0
%                 ofltime = ofltime + WRAPAROUND;
%                 cnt_Ofl = cnt_Ofl+1;
%                 %f1.write("Ofl "+" ")
%             else
%                 cnt_M=cnt_M+1;
%                 %f1.write("MA:%1u "+markers+" ")
%             end
%         end
%             
%         truensync = ofltime + nsync;
%         truetime = (truensync * syncperiod) + (dtime*Resolution);
%         trueTimeArr(b) = truetime;
%         dTimeArr(b) = dtime;
%     end

if strcmp(processFlag, 'Parallel');

    markers = ones(Records, 1);
        
    % Parallelized version
    % Details below, but 1e5 or 1e6 loops and repeated calls is VERY slow. 
    % Can do in parallel here and increase speed
    % Increase in speed 50x for ~1e6 records
    AllRecords = fread(f, Records, 'ubit32');
    nSync = bitand(AllRecords, 65535);
    Channel = bitand(bitshift(AllRecords, -28), 15);
    dTime = bitand(bitshift(AllRecords,-16),4095);
    
    markers(Channel == 15) = bitand(bitshift(AllRecords(Channel == 15),-16),15); %
    overflows = markers == 0;
    ofltime = cumsum(overflows)*WRAPAROUND;
    
    truensync = ofltime + nSync;
    AbsoluteTime = truensync*syncperiod + dTime*Resolution;
    
elseif strcmp(processFlag, 'PQ')
    
    nSync = zeros(Records, 1);
    Channel = zeros(Records, 1);
    dTime = zeros(Records, 1);
    AbsoluteTime = zeros(Records, 1);
    
    
    for k = 1:Records

        % https://github.com/PicoQuant/PicoQuant-Time-Tagged-File-Format-Demos/blob/master/legacy/PicoHarp%20filedemos/pt3/matlab/read_pt3.m
        
        T3Record = fread(f, 1, 'ubit32');     % all 32 bits:
  
    %   +-------------------------------+  +-------------------------------+ 
    %   |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %   +-------------------------------+  +-------------------------------+  

        nSync(k) = bitand(T3Record,65535);       % the lowest 16 bits:
  
    %   +-------------------------------+  +-------------------------------+ 
    %   | | | | | | | | | | | | | | | | |  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
    %   +-------------------------------+  +-------------------------------+  
  
        Channel(k) = bitand(bitshift(T3Record,-28),15);   % the upper 4 bits:

    %   +-------------------------------+  +-------------------------------+ 
    %   |x|x|x|x| | | | | | | | | | | | |  | | | | | | | | | | | | | | | | |
    %   +-------------------------------+  +-------------------------------+


%         dtime=0;
        
        if (Channel(k) < 7) && (Channel(k) > 0)
    
        % if chan = 1 through 6, then these  bits contain the dtime:
        %   +-------------------------------+  +-------------------------------+ 
        %   | | | | |x|x|x|x|x|x|x|x|x|x|x|x|  | | | | | | | | | | | | | | | | |
        %   +-------------------------------+  +-------------------------------+    
        
            dTime(k) = bitand(bitshift(T3Record,-16),4095); 
            
        elseif (Channel(k) == 15)
            
                                                            % This means we have a special record
             markers = bitand(bitshift(T3Record,-16),15); % where these four bits are markers:
     
            %   +-------------------------------+  +-------------------------------+ 
            %   | | | | | | | | | | | | |x|x|x|x|  | | | | | | | | | | | | | | | | |
            %   +-------------------------------+  +-------------------------------+
    
            if markers==0                           % then this is an overflow record
                ofltime = ofltime + WRAPAROUND;         % and we unwrap the numsync (=time tag) overflow
                cnt_Ofl=cnt_Ofl+1;
%                 fprintf(fpout,'Ofl ');
            else                                    % if nonzero, then this is a true marker event
%                 fprintf(fpout,'MA:%1u', markers);
                cnt_M=cnt_M+1;
            end

        end

        truensync = ofltime + nSync(k);
        AbsoluteTime(k) = truensync*syncperiod + dTime(k)*Resolution;
     
    end
    
else
    fprintf(1, 'Optional second argument must be "Parallel" or "PQ"\n');
    
end


    MacroTime = nSync;
    MicroTime = dTime;
    
    fclose(f);
    %f1.close();