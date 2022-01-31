function OUT = morphoTide(input,lat,varargin)
% Morphological tide toolbox (morphoTide)
% Tool box returns a representative tidal signal determined from a tidal
% input signal. The input signal needs to be a scalar vector. This can be,
% for example, tidal water levels or velocity magnitudes.
% Function requires equidistant time-interval of the input signal.
%
% Options for the output are a representative spring-neap cycle
% (Schrijvershof et al, 2022) or a representative double tide (Lesser, 2009)
%
%
% INPUT:
% input     Mx2 numeric array [time,val] consisting of
%           - time: array of matlab datenum values
%           - val:  array of values (water level heights or velocity magnitudes)
%           OR
%           Mx3 cell array consisting of
%           - Constituents names (e.g. M2, M4 etc.)
%           - Constituents amplitudes (in m)
%           - Constituents phases (in degree)
% lat       Latitude (for nodal corrections in tidal analysis)
% Optional:
% type      'springneap' or 'doubletide'
% typeSpec  '' or 'const'
% tStart    Moment to start the representative signal (in datenum)
% timeSpan  Timespan to be modelled in  morphological years
%           (e.g. 10 years)
% nCyles    Number of cycles (repetitions) of the representative signal
% morfac    Morphological acceleration factor used
% plot      0 or 1 (default 0)
%
%
% OUTPUT:
% OUT       representative tidal signal
%
%
%
% References
%
% Lesser, G. R. (2009).
%   An Approach to Medium-term Coastal Morphological Modelling. Retrieved from http://www.narcis.nl/publication/RecordID/oai:tudelft.nl:uuid:62caa573-4fc0-428e-8768-0aa47ab612a9
% Roelvink, D., & Reniers, A. (2011).
%   A Guide to Modeling Coastal Morphology. WORLD SCIENTIFIC. https://doi.org/10.1142/7712


%% Settings

OPT.type        = 'springneap';
OPT.typeSpec    = '';
OPT.tStart      = datenum(2017,12,25);
OPT.tStop       = datenum(2019,01,01);
OPT.timeSpan    = [];
OPT.nCycles     = 2;
OPT.morfac      = [];
OPT.plot        = 0;

if nargin > 2
    for i = 1:length(varargin)
        if isfield(OPT,varargin{i})
            fld = varargin{i};
            OPT.(fld) = varargin{i+1};
        end
    end
end

%% Define and rewrite variables

% Set-up struct
OUT = struct;

%% Check the input

% Time vector (days, hours, and seconds)
if isnumeric(input)
    OUT.td = (input(1,1):datenum(0,0,0,0,10,0):input(end,1))';
elseif iscell(input)
    OUT.td = (OPT.tStart:datenum(0,0,0,0,10,0):OPT.tStart+OPT.nCycles*14.77)';
end
OUT.th      = OUT.td*24;
OUT.ts      = OUT.td*24*60*60;


% Tidal analysis and rewrite tidal constituents
const = [];
if isnumeric(input)
    
    OUT.valOri = input(:,2);
    
    % Tidal analysis (requires t_tide toolbox)
    if ~exist('t_tide.m','file')
        fprintf('\tThe morphoTide toolbox requires the t_tide toolbox version 1.4,\n\tdownload it at https://www.eoas.ubc.ca/~rich/\n');
    end
    
    
    T = struct;
    try
        % Harmonic analysis analysis (t_tide)
        % Version 1.4: Use constant value for time interval
        [OUT.T,OUT.valTide] = t_tide(input(:,2),...
            'int',diff(input(1:2,1)).*24,...
            'start',input(1,1),...
            'latitude',lat,...
            'output','none');
        OUT.T.names = cellstr(OUT.T.name);
        OUT.T.per   = 1./OUT.T.freq;  % period (hr)
        OUT.T.cel   = OUT.T.freq.*360; % celerity (deg/hr)
        % Significant peaks
        OUT.T.fsig = OUT.T.tidecon(:,1) > OUT.T.tidecon(:,2);
    catch
        fprintf('\tAn old version of t_tide is present,\n\tplease download the most recent version (1.4) at https://www.eoas.ubc.ca/~rich/,\n\tand run morphoTide again.\n');
        return
    end
    
    const.name = cellstr(OUT.T.name(OUT.T.fsig,:));
    const.amp  = OUT.T.tidecon(OUT.T.fsig,1);
    const.pha  = OUT.T.tidecon(OUT.T.fsig,3);
    const.freq = OUT.T.freq(OUT.T.fsig); % cycles/hr
    
elseif iscell(input)
    const.name  = input(:,1);
    const.amp   = cell2mat(input(:,2));
    const.pha   = cell2mat(input(:,3));
    const.freq  = t_tide_name2freq(const.name,'unit','cyc/hr')';
end


%% Representative signal

switch OPT.type
    case 'springneap'
        
        
        % Reconstruct a time-series from the given tidal components
        if isnumeric(input)
            OUT.valOri = input(:,2);
        elseif iscell(input)
            tidecon = eps*ones(length(const.name),4);
            tidecon(:,1) = const.amp;
            tidecon(:,3) = const.pha;
            yout = t_predic(OUT.td,const.name,const.freq,tidecon,...
                'start time',[]); % Without nodal corrections and no adjustment to central time
            % rewrite
            OUT.valOri      = yout(:);
            OUT.valTide     = yout(:);
        end
        
        
        if isempty(OPT.typeSpec)
            
            %% Power spectrum
            % Remove NaN
            nnid    = isnan(OUT.valOri);
            val2    = interp1(OUT.td(~nnid),OUT.valOri(~nnid),OUT.td);
            dt      = round(diff(OUT.td(1:2))*24*60*60);  % Sampling period (s)
            Fs      = 1/dt;                      % sampling frequency (Hz)
            
            nfft = length(val2);
            fftx = fft(val2,nfft);
            P2   = abs(fftx/nfft);
            L    = floor(nfft/2+1);
            P1   = P2(1:L,1);
            
            fHz = Fs*(0:(nfft/2))/nfft; fHz = fHz(:);
            fHr = fHz/(1/3600); % Frequency for tides more common in hr^-1
            
            
            % 1. Find amplitudes of diurnal (D1), semi-diurnal (D2) and quarter diurnal
            %    (D4) frequencies
            %    D1: periods (hr) 21.5782 (UPS1) - 29.0727 (ALP1)
            %    D2: periods (hr) 11.4067 (2SN2) - 13.6323 (ST36)
            %    D3: periods (hr)  7.9927 (SK2)  - 8.4940 (NO3)
            %    D4: periods (hr)  5.9918 (SK4)  - 6.3821 (ST8)
            %    D6: periods (hr)  4.2154 (ST11) - 4 (S6)
            %    D8: periods (hr)  3.076  (3MK8) - 3.1346 (ST18)
            
            %         name        deg/hr        hr/cyc        cyc/hr       cyc/day
            %         O1          13.943       25.8193      0.038731       0.92954
            %         K1         15.0411       23.9345      0.041781        1.0027
            
            SN = struct;
            % Define frequency bands for harmonic analysis
            SN.D1.fmin      = 1/29.0727; % Diurnals
            SN.D1.fmax      = 1/21.5782;
            SN.D2.fmin      = 1/13.6323; % Semi-diurnals
            SN.D2.fmax      = 1/11.4067;
            SN.D4.fmin      = 1/ 6.3821; % Quarter-diurnals
            SN.D4.fmax      = 1/ 5.9918;
            SN.D6.fmin      = 1/4.2154;
            SN.D6.fmax      = 1/4;
            SN.D8.fmin      = 1/3.1346;
            SN.D8.fmax      = 1/3.076;
            
            
            % Find peak frequencies in the defined bands
            flds = {'D1','D2','D4','D6','D8'};
            for i = 1:length(flds)
                f = flds{i}; clear pwr
                idfb = fHr > SN.(f).fmin & fHr < SN.(f).fmax;
                pwr         = P1;
                pwr(~idfb)  = 0;
                
                [pks,idpks] = findpeaks(pwr,fHr,'MinPeakDistance',0.001);
                
                [pks,ids] = sort(pks,'descend');
                idpks     = idpks(ids);
                
                % First peak
                % D1: O1 or K1
                % D2: M2
                % D4: M4
                SN.(f).fpeak1    = idpks(1);
                SN.(f).Tpeak1    = 1/SN.(f).fpeak1;
                SN.(f).pwrpeak1  = pks(1);
                % Second peak
                % D1: O1 or K1
                % D2: S2 or N2
                % D4: M6
                if length(ids) > 1
                    SN.(f).fpeak2    = idpks(2);
                    SN.(f).Tpeak2    = 1/SN.(f).fpeak2;
                    SN.(f).pwrpeak2  = pks(2);
                end
            end
            
            tc = t_getconsts;
            % The found peak frequencies sometimes show sligth differences
            % from the main tidal frequencies expected, perhaps because of
            % the 'MinPeakDistance', solved now by setting the peak
            % frequencies to tidal frequencies if the difference is less
            % than 'MinPeakDistance'.
            for i = 1:length(flds)
                f = flds{i};
                [min1,id1] = min(abs(tc.freq-SN.(f).fpeak1));
                if min1 < 0.001
                    SN.(f).fpeak1 = tc.freq(id1);
                    SN.(f).Tpeak1 = 1/SN.(f).fpeak1;
                end
                
                if isfield(SN.(f),'fpeak2')
                    [min2,id2] = min(abs(tc.freq-SN.(f).fpeak2));
                    if min2 < 0.001
                        SN.(f).fpeak2 = tc.freq(id2);
                        SN.(f).Tpeak2 = 1/SN.(f).fpeak2;
                    end
                end
            end
            
            
            
            
            %% Harmonic analysis on the peak frequencies
            
            periods(1,1)    = SN.D1.Tpeak1*60*60;
            periods(2,1)    = SN.D1.Tpeak2*60*60;
            periods(3,1)    = SN.D2.Tpeak1*60*60;
            periods(4,1)    = SN.D2.Tpeak2*60*60;
            periods(5,1)    = SN.D4.Tpeak1*60*60;
            periods(6,1)    = SN.D4.Tpeak2*60*60;
            periods(7,1)    = SN.D6.Tpeak1*60*60;
            periods(8,1)    = SN.D6.Tpeak2*60*60;
            periods(9,1)    = SN.D8.Tpeak1*60*60;
            
            FIT = harmanal(OUT.ts,OUT.valOri,'T',periods);
            
            % Rewrite to constituents
            n = 0;
            for i = 1:length(flds)
                f = flds{i};
                for j = 1:2
                    n = n+1;
                    if n > length(FIT.hamplitudes)
                        continue;
                    end
                    SN.(f).A(j,1)        = FIT.hamplitudes(n);
                    SN.(f).phi(j,1)      = FIT.hphases(n);
                    SN.(f).omega(j,1)    = FIT.omega(n);
                    SN.(f).T(j,1)        = FIT.T(n);
                end
            end
            
        elseif strcmp(OPT.typeSpec,'const')
            % Perform the same analysis directly from the tidal
            % constituents
            
            idO1 = strcmp(const.name,'O1');
            idK1 = strcmp(const.name,'K1');
            idM2 = strcmp(const.name,'M2');
            idS2 = strcmp(const.name,'S2');
            idM4 = strcmp(const.name,'M4');
            idM6 = strcmp(const.name,'M6');
            idM8 = strcmp(const.name,'M8');
            
            SN.D1.A(1,1) = const.amp(idO1);
            SN.D1.A(2,1) = const.amp(idK1);
            SN.D2.A(1,1) = const.amp(idM2);
            SN.D2.A(2,1) = const.amp(idS2);
            SN.D4.A(1,1) = const.amp(idM4);
            SN.D6.A(1,1) = const.amp(idM6);
            SN.D8.A(1,1) = const.amp(idM8);
            
            
            SN.D1.phi(1,1) = deg2rad(const.pha(idO1));
            SN.D1.phi(2,1) = deg2rad(const.pha(idK1));
            SN.D2.phi(1,1) = deg2rad(const.pha(idM2));
            SN.D2.phi(2,1) = deg2rad(const.pha(idS2));
            SN.D4.phi(1,1) = deg2rad(const.pha(idM4));
            SN.D6.phi(1,1) = deg2rad(const.pha(idM6));
            SN.D8.phi(1,1) = deg2rad(const.pha(idM8));
            
            SN.D1.T(1,1) = (1/const.freq(idO1))*60*60;
            SN.D1.T(2,1) = (1/const.freq(idK1))*60*60;
            SN.D2.T(1,1) = (1/const.freq(idM2))*60*60;
            SN.D2.T(2,1) = (1/const.freq(idS2))*60*60;
            SN.D4.T(1,1) = (1/const.freq(idM4))*60*60;
            SN.D6.T(1,1) = (1/const.freq(idM6))*60*60;
            SN.D8.T(1,1) = (1/const.freq(idM8))*60*60;
            
            SN.D1.omega(1,1) = (2*pi)/SN.D1.T(1,1);
            SN.D1.omega(2,1) = (2*pi)/SN.D1.T(2,1);
            SN.D2.omega(1,1) = (2*pi)/SN.D2.T(1,1);
            SN.D2.omega(2,1) = (2*pi)/SN.D2.T(2,1);
            SN.D4.omega(1,1) = (2*pi)/SN.D4.T(1,1);
            SN.D6.omega(1,1) = (2*pi)/SN.D6.T(1,1);
            SN.D8.omega(1,1) = (2*pi)/SN.D8.T(1,1);
        else
            fprintf('\tVariable typeSpec not recognized\n');
        end
        
        
        
        
        %% Build-up of representative periodic spring-neap cycle
        
        
        % 1. Determine duration of the synthetic signal (spring-neap cyle)
        % 2. Limit the length to an integer and even number of D2 cycles
        % 3. Create a time vector with second increments
        
        % 1.
        dur_s = ((2*pi) / (SN.D2.omega(2) - SN.D2.omega(1)));
        dur_d = dur_s/60/60/24;
        
        % 2.
        nD2cycles   = floor(dur_s/SN.D2.T(1)/2)*2;
        Tmor        = nD2cycles*SN.D2.T(1);
        
        % 3.
        SN.mor.ts = (0:Tmor)';
        SN.mor.td = SN.mor.ts/60/60/24;
        
        % Caclulate a component Dsn that mimics the spring-neap variation
        % during the cyle
        SN.Dsn.omega    = (2*pi)/Tmor;
        SN.Dsn.Avar     = SN.D2.A(2) * cos(SN.Dsn.omega.*SN.mor.ts - pi); % -pi to start with neap tide
        SN.Dsn.Avar2    = SN.D2.A(2) * cos(SN.Dsn.omega.*OUT.ts - SN.D2.phi(2)-pi); % to synchronize with original signal
        
        % To check: Very small rounding differences can cause the end of
        % the singal not to coicide with the beginning of the signal. This
        % is solved by making Avar exactly symmetric
        SN.Dsn.Avar  = round(SN.Dsn.Avar,4);
        
        
        % Calculate C1 component from O1 & K1 (amplitude and phase)
        % Following Lesser et al. (2009)
        D1acos = cos(SN.D1.phi(1));
        D1asin = sin(SN.D1.phi(1));
        D1bcos = cos(SN.D1.phi(2));
        D1bsin = sin(SN.D1.phi(2));
        D1aim  = D1acos + 1i*D1asin;
        D1bim  = D1bcos + 1i*D1bsin;
        pha    = (D1aim+D1bim)/2;
        
        SN.C1.omega = 0.5*SN.D2.omega(1);
        SN.C1.A     = sqrt(2*SN.D1.A(1)*SN.D1.A(2));
        SN.C1.phi   = mod(atan2d(imag(pha),real(pha)),360);   % to check
        %         SN.C1.A = SN.D1.A(1);
        %         SN.C1.phi = SN.D1.phi(1);
        
        % Amplification of M2, following Lesser (2009) equation 5.16
        % sqrt((Aconst^2 - Ac1^2) / Am2^2)
        SN.D2.ampfac = sqrt((sum(const.amp.^2)-SN.C1.A^2)/SN.D2.A(1)^2);
        SN.D2.A(1) = SN.D2.A(1)*SN.D2.ampfac;
        
        
        
        %%% Step-wise set-up of the representative spring-neap cycle %%%
        
        % 1. Diurnal, semi-diurnal, and quarter-diurnal signal
        SN.mor.hD1a = ...
            (SN.D1.A(1)) .* cos(SN.D1.omega(1).*SN.mor.ts - SN.D1.phi(1));
        SN.mor.hD1b = ...
            (SN.D1.A(2)) .* cos(SN.D1.omega(2).*SN.mor.ts - SN.D1.phi(2));
        SN.mor.hC1 = ...
            (SN.C1.A) .* cos(SN.C1.omega.*SN.mor.ts - SN.C1.phi);
        SN.mor.hD2 = ...
            (SN.D2.A(1)) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1));
        SN.mor.hD4 = ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1));
        SN.mor.hD4D6D8 = ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
            (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
        
        % 2. Inc. mimic of spring-neap variation
        SN.mor.hD2mod = ...
            (SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1));
        % 2b Spring-neap variation through M2S2 interaction
        SN.mor.hM2S2 = ...
            (SN.D2.A(1) .* cos(SN.D2.omega(1).*SN.mor.ts)) + ...
            (SN.D2.A(2) .* cos(SN.D2.omega(2).*SN.mor.ts));
        % Including quarter-diurnals
        SN.mor.hD2modD4 = ...
            (SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1));
        % Including 6th diurnals
        SN.mor.hD2modD4D6 = ...
            (SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1));
        % Including 8th diurnals
        SN.mor.hD2modD4D6D8 = ...
            (SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
            (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
        % Including diurnals component D1
        SN.mor.hD1D2modD4D6D8 = ...
            (SN.D1.A(1)) .* cos(SN.D1.omega(1).*SN.mor.ts - SN.D1.phi(1)) + ...
            (SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
            (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
        % Including diurnals component C1
        SN.mor.hC1D2modD4D6D8 = ...
            (SN.C1.A) .* cos(SN.C1.omega.*SN.mor.ts - SN.C1.phi) + ...
            (SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
            (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
        
        %% Scaling
        % Scaling of semi-diurnal amplitude
        val1 = OUT.valTide;
        dx        = 0.2; % Bin interval
        fac = 0.5:0.01:1.5;
        ERR = NaN(length(fac),1);
        for k = 1:length(fac)
            val2 = ...
                (SN.C1.A) .* cos(SN.C1.omega.*SN.mor.ts - SN.C1.phi) + ...
                (fac(k) * SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
                (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
                (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
                (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
            edges   = round(min([val1;val2]),1,'decimals')-dx:dx:round(max([val1;val2]),1,'decimals')+dx;
            bins    = edges(1:end-1)+diff(edges);
            data1 = histcounts(val1,edges,'Normalization','probability');
            data2 = histcounts(val2,edges,'Normalization','probability');
            % RMSE
            ERR(k,1) = sqrt( sum((data2-data1).^2) / length(data1));
        end
        [~,id] = min(ERR);
        fD2 = fac(id);
        SN.mor.hC1fD2modD4D6D8 = ...
            (SN.C1.A) .* cos(SN.C1.omega.*SN.mor.ts - SN.C1.phi) + ...
            (fD2 * SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
            (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
        % Scaling of quarter-diurnal amplitude
        ERR = NaN(length(fac),1);
        for k = 1:length(fac)
            val2 = ...
                (SN.C1.A) .* cos(SN.C1.omega.*SN.mor.ts - SN.C1.phi) + ...
                (fD2 * SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
                (fac(k) * SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
                (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
                (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
            
            edges   = round(min([val1;val2]),1,'decimals')-dx:dx:round(max([val1;val2]),1,'decimals')+dx;
            bins    = edges(1:end-1)+diff(edges);
            data1 = histcounts(val1,edges,'Normalization','probability');
            data2 = histcounts(val2,edges,'Normalization','probability');
            % RMSE
            ERR(k,1) = sqrt( sum((data2-data1).^2) / length(data1));
        end
        [~,id] = min(ERR);
        fD4 = fac(id);
        SN.mor.hC1fD2modfD4D6D8 = ...
            (SN.C1.A) .* cos(SN.C1.omega.*SN.mor.ts - SN.C1.phi) + ...
            (fD2 * SN.D2.A(1) + SN.Dsn.Avar) .* cos(SN.D2.omega(1).*SN.mor.ts - SN.D2.phi(1)) + ...
            (fD4 * SN.D4.A(1)) .* cos(SN.D4.omega(1).*SN.mor.ts - SN.D4.phi(1)) + ...
            (SN.D6.A(1)) .* cos(SN.D6.omega(1).*SN.mor.ts - SN.D6.phi(1)) + ...
            (SN.D8.A(1)) .* cos(SN.D8.omega(1).*SN.mor.ts - SN.D8.phi(1));
        
        %% Create a signal with n cycles of repetition
        % 1. Concatenate the constructed signal
        % 2. Resample to 10 min intervals
        
        clear val;
        val         = repmat(SN.mor.hC1fD2modfD4D6D8,OPT.nCycles,1);
        ts          = (1:1:length(val))';
        SN.mor.val  = val(1:600:end);
        ts          = ts(1:600:end);
        td          = ts/60/60/24;
        
        
        OUT.SN      = SN;
        OUT.output  = SN.mor.val;
        OUT.datenum = OPT.tStart+td;
        
        
        tFull = OUT.td(end) - OUT.td(1);
        tSN   = OUT.SN.mor.td(end) - OUT.SN.mor.td(1);
        n     = floor(tFull/tSN);
        
        val     = repmat(SN.mor.hC1fD2modfD4D6D8,n,1);
        OUT.val = NaN(length(OUT.td),1);
        tmp     = val(1:600:end);
        OUT.val(1:length(tmp)) = tmp;
        
    case 'doubletide'
        
        
        clear D
        % Constituents for analysis and rewrite significant constituents to struct
        D.const  = {'M2','M4','M6','M8','M10','O1','K1'};
        D.names  = const.name;
        D.amps   = const.amp;
        D.phas   = const.pha;
        D.phas   = mod(D.phas,360.);
        D.freqs  = const.freq; % cyc/hr
        D.cels   = D.freqs.*360; % celerity (deg/hr)       
        
        % find constituent and rewrite to separate fields
        for i = 1:length(D.const)
            c = D.const{i};
            if any(strcmp(D.names,c))
                id = strcmp(D.names,c);
                D.mor.(c).amp     = D.amps(id);
                D.mor.(c).pha     = D.phas(id);
                D.mor.(c).freq    = D.freqs(id);
                D.mor.(c).cel     = D.cels(id);
                D.mor.(c).vel     = 1/D.mor.(c).cel*360*60;  % Velocity
            end
        end
        
        % Calculate C1 component from O1 & K1 (amplitude and phase)
        if ~isfield(D.mor,'O1')
            fprintf('\tO1 constituent could not be determined!\n');
        elseif ~isfield(D.mor,'K1')
            fprintf('\tK1 constituent could not be determined!\n');
        end
        D.mor.C1.amp    = sqrt(2*D.mor.O1.amp*D.mor.K1.amp);
        D1acos          = cosd(D.mor.O1.pha);
        D1asin          = sind(D.mor.O1.pha);
        D1bcos          = cosd(D.mor.K1.pha);
        D1bsin          = sind(D.mor.K1.pha);
        D1aim           = D1acos + 1i*D1asin;
        D1bim           = D1bcos + 1i*D1bsin;
        pha             = (D1aim+D1bim)/2;
        
        D.mor.C1.pha    = mod(atan2d(imag(pha),real(pha)),360);   % to check
        D.mor.C1.freq   = 0.5*D.mor.M2.freq;
        D.mor.C1.cel    = 0.5*D.mor.M2.cel;
        D.mor.C1.vel    = 1/D.mor.C1.cel*360*60;
        
        % Amplification of M2, following Lesser (2009) equation 5.16
        % sqrt((Aconst^2 - Ac1^2) / Am2^2)
        D.ampfac = sqrt((sum(D.amps.^2)-D.mor.C1.amp.^2)/D.mor.M2.amp.^2);
        D.mor.M2.amp = D.mor.M2.amp*D.ampfac;
        
        
        % Rewrite
        D.constMor = fieldnames(D.mor);
        D.mor.const = {'M2','M4','M6','M8','M10','C1'};
        n = 0;
        for i = 1:length(D.mor.const)
            c = D.mor.const{i};
            if any(strcmp(D.constMor,c))
                n = n+1;
                D.mor.amp(n,1) = D.mor.(c).amp;
                D.mor.pha(n,1) = D.mor.(c).pha;
                D.mor.freq(n,1) = D.mor.(c).freq;
                D.mor.cel(n,1) = D.mor.(c).cel;
                D.mor.vel(n,1) = D.mor.(c).vel;
            end
        end
        D.mor.per_h = 1./D.mor.freq;
        D.mor.per_s = D.mor.per_h*3600;
        D.mor.omega = (2*pi)./(D.mor.per_h); % deg/hr

        for i = 1:length(D.mor.amp)
            pha = deg2rad(D.mor.pha(i));
            D.mor.val(:,i) = D.mor.amp(i) * cos(D.mor.omega(i)*OUT.th + pha);
        end
        D.mor.eta = sum(D.mor.val,2);
        
        % Write to output struct
        OUT.val     = D.mor.eta(:);
        OUT.output  = [D.mor.vel,D.mor.amp,D.mor.pha];
        
    otherwise
        fprintf('\tNo type defined\n');
end

%% Figures settings




%% Figure of power spectrum
if OPT.plot
    
% Some general plotsettings
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex');

    
    %                         close all;
    %                         fig = figure; fig.Units = 'centimeters'; fig.Position = [5 5 12 10];
    %                         axs = tight_subplot(1,1,0.5,[0.15,0.05],[0.15,0.05]);
    %                         set(fig,'CurrentAxes',axs(1)); ax = gca; hold on; box on;
    %
    %                         plot(fHr,P1,'.-k');
    %                         for i = 1:length(flds)
    %                             f = flds{i};
    %                             text(SN.(f).fpeak1,SN.(f).pwrpeak1,f,...
    %                                 'HorizontalAlignment','c','VerticalAlignment','bottom')
    %                         end
    %                         xlim([0,1]); ylim([1e-8,10])
    %                         ax.XScale = 'log';
    %                         ax.YScale = 'log';
    %                         xlabel('$f (hr^{-1})$');
    %                         ylabel('Power');
    
    %% Figure of time-series and histograms
    clr(1,:)    = [0,0,0];
    clr(2:3,:)  = lines(2);
    
    % Rewrite
    if isnumeric(input)
        val0 = OUT.valOri;
    else
        val0 = NaN;
    end
    val1 = OUT.valTide;
    val2 = OUT.val;
%     val3 = OUT.SN.mor.hC1fD2modfD4D6D8;
    
    
    close all; clear axs
    % Set-up figure
    fig = figure; fig.Units = 'centimeters'; fig.Position = [-40 0 34 19];
    axs(1) = axes; axs(1).Units = 'centimeters'; axs(1).Position = [1.5 10 22 8.5];
    axs(2) = axes; axs(2).Units = 'centimeters'; axs(2).Position = [25 10 8 8.5];
    axs(3) = axes; axs(3).Units = 'centimeters'; axs(3).Position = [1.5 1.5 9 7];
    axs(4) = axes; axs(4).Units = 'centimeters'; axs(4).Position = [12.5 1.5 5 7];
    axs(5) = axes; axs(5).Units = 'centimeters'; axs(5).Position = [18.0 1.5 5 7];
    axs(6) = axes; axs(6).Units = 'centimeters'; axs(6).Position = [24.5 1.5 9 7];
    
    
    set(fig,'CurrentAxes',axs(1)); hold on; box on; ax = gca;
    plot(OUT.td,val0,'--','Color',clr(1,:));
    plot(OUT.td,val1,'Color',clr(1,:));
    plot(OUT.td,val2,'Color',clr(2,:));
    
    axis tight;
    xlim([OUT.td(1),OUT.td(end)])
    datetick('x','keeplimits')
    legStr = {'Measured','Full tidal','Morph. spring-neap (non-scaled)','Morph. spring-neap (scaled)'};
    hl = legend(legStr,'Location','NE');
    
    set(fig,'CurrentAxes',axs(2)); hold on; box on; ax = gca;
    plot(OUT.td,val0,'--','Color',clr(1,:));
    plot(OUT.td,val1,'Color',clr(1,:));
    plot(OUT.td,val2,'Color',clr(2,:));
    
    axis tight;
    xlim([OUT.td(1),OUT.td(1)+1])
    datetick('x','HH:mm','keeplimits')
    
    set(fig,'CurrentAxes',axs(3)); hold on; box on; ax = gca;
    dx      = 0.2; % Bin interval
    edges   = round(min([val0;val1;val2]),1,'decimals')-2*dx:dx:round(max([val0;val1;val2]),1,'decimals')+2*dx;
    bins    = edges(1:end-1)+diff(edges);
    
    data1 = histcounts(val1,edges,'Normalization','probability');
    data2 = histcounts(val2,edges,'Normalization','probability');
    
    area(bins,data1,'LineStyle','-','EdgeColor','n',...
        'LineWidth',1,'FaceColor',clr(1,:),'FaceAlpha',0.5);
    plot(bins,data2,'Color',clr(2,:));
    
    xlabel('\zeta (m)')
    
    DUR = tideTools.duration(OUT.td,val1);
    val1rt = DUR.rt;
    val1ft = DUR.ft;
    DUR = tideTools.duration(OUT.td,val2);
    val2rt = DUR.rt;
    val2ft = DUR.ft;

    
    set(fig,'CurrentAxes',axs(4)); hold on; box on; ax = gca;
    % Settings
    dx      = 0.25; % Bin interval
    edges   = 4:dx:9;
    bins    = edges(1:end-1)+diff(edges);
    norm    = 'probability';
    
    data1 = histcounts(val1ft.*24,edges,'Normalization',norm);
    data2 = histcounts(val2ft.*24,edges,'Normalization',norm);

    area(bins,data1,'LineStyle','-','EdgeColor','n',...
        'LineWidth',1,'FaceColor',clr(1,:),'FaceAlpha',0.5);
    plot(bins,data2,'Color',clr(2,:));
    xlim([min(edges),max(edges)]);
    xlabel('Dur. falling tide (hrs)')
    
    set(fig,'CurrentAxes',axs(5)); hold on; box on; ax = gca;
    data1 = histcounts(val1rt.*24,edges,'Normalization',norm);
    data2 = histcounts(val2rt.*24,edges,'Normalization',norm);
    
    area(bins,data1,'LineStyle','-','EdgeColor','n',...
        'LineWidth',1,'FaceColor',clr(1,:),'FaceAlpha',0.5);
    plot(bins,data2,'Color',clr(2,:));
    xlim([min(edges),max(edges)]);
    xlabel('Dur. rising tide (hrs)')

    set(fig,'CurrentAxes',axs(6)); hold on; box on; ax = gca;
    dzdt1 = diff(val1)./(diff(OUT.td)*24);
    dzdt2 = diff(val2)./(diff(OUT.td)*24);
    
    dx      = 1/6; % Bin interval
    edges   = round(min([dzdt1;dzdt2]),1,'decimals')-2*dx:dx:round(max([dzdt1;dzdt2]),1,'decimals')+2*dx;
    bins    = edges(1:end-1)+diff(edges);
    
    data1 = histcounts(dzdt1,edges,'Normalization','probability');
    data2 = histcounts(dzdt2,edges,'Normalization','probability');
    
    area(bins,data1,'LineStyle','-','EdgeColor','n',...
        'LineWidth',1,'FaceColor',clr(1,:),'FaceAlpha',0.5);
    plot(bins,data2,'Color',clr(2,:));
    
    xlabel('$d\zeta/dt (m/hr)$');
    
    
end

return




