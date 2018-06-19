% Generate area-normalized reference spectra from .pt3 files

LoFilepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\Figure1\LUVsB_PSM_Chol_NR12S_2.pt3';
LdFilepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\Figure1\SLBs_DOPC_Chol_NR12S_3.pt3';


outputFilepath = 'D:\Dropbox\Proposals\FSCS\DataFiles\ReferenceSpectra\';

% Channel ranges pulled from .fcs file metadata and summarized here
ChannelRanges = [535.03 570.67;
                 570.67 588.48;
                 588.48 606.30;
                 606.30 624.11;
                 624.11 650.84;
                 650.84 695.38];
             
% Channel centers
ChanCents = mean(ChannelRanges, 2);

% Lo phase
[Chan, ~, ~, MicroTime] =  pt3Import(LoFilepath);

LoSpectrum = histc(MicroTime(Chan < 15), 0:5);
LoSpectrum = LoSpectrum/sum(LoSpectrum);

% Ld phase
[Chan, ~, ~, MicroTime] =  pt3Import(LdFilepath);

LdSpectrum = histc(MicroTime(Chan < 15), 0:5);
LdSpectrum = LdSpectrum/sum(LdSpectrum);

% Save reference spectra to ASCII files to be used in other measurements
[~, Lofn] = fileparts(LoFilepath);
[~, Ldfn] = fileparts(LdFilepath);

LoFID = fopen(fullfile(outputFilepath, strcat(Lofn, '_spectrum.txt')), 'w+');
fprintf(LoFID, '# FSCS Spectrum File - created by GenerateReferenceSpectra.m\r\n');
fprintf(LoFID, '# Data File - %s\r\n', LoFilepath);
fprintf(LoFID, '# # # # # # # # # # # # # # # # # # # # # # # #\r\n');
fprintf(LoFID, 'ChannelCenter(nm)\tIntensity(au)\r\n');
for k = 1:numel(ChanCents);
    fprintf(LoFID, '%.4f\t%.4f\r\n', [ChanCents(k) LoSpectrum(k)]);
end
fclose(LoFID);

LdFID = fopen(fullfile(outputFilepath, strcat(Ldfn, '_spectrum.txt')), 'w+');
fprintf(LdFID, '# FSCS Spectrum File - created by GenerateReferenceSpectra.m\r\n');
fprintf(LdFID, '# Data File - %s\r\n', LoFilepath);
fprintf(LdFID, '# # # # # # # # # # # # # # # # # # # # # # # #\r\n');
fprintf(LdFID, 'ChannelCenter(nm)\tIntensity(au)\r\n');
for k = 1:numel(ChanCents);
    fprintf(LdFID, '%.4f\t%.4f\r\n', [ChanCents(k) LdSpectrum(k)]);
end
fclose(LdFID);