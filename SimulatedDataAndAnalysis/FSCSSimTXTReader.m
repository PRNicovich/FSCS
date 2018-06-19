% Reader for FSCS simulation output text files

function dataReturn = FSCSSimTXTReader(fileName)

    fID = fopen(fileName, 'r');
    eofCheck = false;

    dataOut = [];
    dataChunk = [];

    while ~eofCheck

        nextLine = fgetl(fID);

        if nextLine == -1
            eofCheck = true;
            dataOut = [dataOut dataChunk];
            break

        else

            if strcmp(nextLine(1:2), '# ') 
                % Headerline
            elseif strcmp(nextLine(1:2), '##')
                % Separator lines
                dataOut = [dataOut dataChunk];
                dataChunk = [];

            else
                % Actual data or column header
                if isstrprop(nextLine(1), 'alpha')
                    % Column header
                else
                    % Actual data line
                    dataLine = textscan(nextLine, '%.1f\t%.6f\t%.6f\t%.6f');
                    dataChunk = [dataChunk; cell2mat(dataLine)];

                end
            end

        end

    end
    
    dataReturn = dataOut;
    
end
        

        