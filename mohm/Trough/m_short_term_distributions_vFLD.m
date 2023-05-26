function[output] = m_short_term_distributions_vFLD(func,data,funcformat)
%
% DESCRIPTION: Collection of short-term distributions conditional on
%              the underlying Hm0. Called by m_mode
% 
% SYNTAX: output = m_short_term_distributions(func,DistIn,flag)
%
% REQUIRED INPUT:   func    string identifier for short-term dist.
%                   DistIn  structure with various input (created in m_mode) 
%                   flag    identifying type of output
%                           1: CDF of individual waves
%                           2: Sample of sea state maximum wave (crest)
%                           3: CDF of maxima
%
% Available short-term distributions are (case sensitive):
%
% Rayleigh_H                Rayleigh wave height distributions
% Forristall_H              Forristall (1978) wave height distributions
% Naess_H                   N?ss (1985) wave height distributions
% Battjes_Groenendijk_H     Battjes & Groenendijk (2000)  wave height distributions
% Glukhovskiy_H             Glukhovskiy wave height distribution, as
%                           presented in Battjes & Groenendijk (2000) paper
% Weibull_H                 A Weibull distribution with user specified
%                           scale (beta) and shape alfa
% Forristall_C              Forristall (2000) crest height distributions
% Rayleigh_T                Rayleigh Trough (added by FLD)
% Arhan_Plaisted_T          Arhan and Plaisted Trough (added by FLD)
%
% Copyright DHI ?, Denmark, All rights reserved
% 2012-06-15 - Hans Fabricius Hansen, hfh@dhigroup.com
%
% Revision History
% HFH 2014-03-20:	Bug fix in B&G distribution of Hmax
% HFH 2018-11-13:   Added funcformat 2 option - sampling of sea state
%                   maximum and some clean up etc.
% FLD 2021-02-03:   Added Rayleigh Trough and Arhan and Plaisted Trough
%                   (use with caution as not extensively tested)
%
func = str2func(func);              % create function handle
output = func(data,funcformat);     % call short-term distribution function
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Rayleigh_H(data,funcformat)
%
% Syntax: CDF = Rayleigh_H(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave heights 
%                 conditional on Hm0
%
% funcformat = 2: Return sample of maximum wave heigth in sea state
%
% funcformat = 3: Return cumulative distribution of Hmax
%                 conditional on Hmp
%
%
%
% perform calculations
switch funcformat
    
    % return cumulative distribution of wave heights conditional on Hm0
    case 1
        H    = reshape(data.X  ,[],1);
        Hm0  = reshape(data.Hm0,1,[]);
        output.CDF = 1 - exp(-2*(H*(Hm0.^-1)).^2);

    % return random sample of max wave in sea state
    case 2
        if ~isfield(data,'P')
            data.P = rand(size(data.Hm0));
        end
        output     = sqrt(2)/2.*data.Hm0.*(-log(1-data.P.^(data.N.^-1))).^(1/2);
        
    % return cumulative distribution of Hmax conditional on Hmp
    case 3
        Hmx  = reshape(data.Xmx,[],1);
        H    = reshape(data.X  ,1,[]);
        N    = data.N;
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^2-1)));
        
    otherwise
        disp('Error in call to Rayleigh distribution function')
        return
end

end

function[output] = Rayleigh_T(data,funcformat)
%
% Syntax: CDF = Rayleigh_T(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave troughs 
%                 conditional on Hm0
%
% based on China Ocean Eng., 2018, Vol. 32, No. 6, P. 655–664 Predicting Nonlinear Wave Trough Distributions Utilizing A Transformed Linear Simulation Method, WANG Ying-guanga https://link.springer.com/article/10.1007/s13344-018-0067-0
%
% perform calculations
switch funcformat
    
    % return cumulative distribution of wave troughs conditional on Hm0
    case 1
        H    = reshape(data.X  ,[],1);
        Hm0  = reshape(data.Hm0,1,[]);
        output.CDF = 1 - exp(-8*(H*(Hm0.^-1)).^2);
    
    % return cumulative distribution of Hmax conditional on Tmp % correct
    % or not??
    case 3
        Hmx  = reshape(data.Xmx,[],1);
        H    = reshape(data.X  ,1,[]);
        N    = data.N;
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^2-1)));
        
    otherwise
        disp('Error in call to Rayleigh T distribution function')
        return
end

end


function[output] = Arhan_Plaisted_T(data,funcformat)
%
% Syntax: CDF = ArhanPlaisted_T(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave troughs 
%                 conditional on Hm0
%
% based on China Ocean Eng., 2018, Vol. 32, No. 6, P. 655–664 Predicting Nonlinear Wave Trough Distributions Utilizing A Transformed Linear Simulation Method, WANG Ying-guanga https://link.springer.com/article/10.1007/s13344-018-0067-0
d     = reshape(data.d  ,1,[]);
T     = reshape(data.Tp,1,[]); % use Tp and T02 as ref?
for tt = 1:length(T)
    lambda(tt) = m_DispRel(d(tt),'T',T(tt)); %SJA took tt out of loop as d is one value
end
% perform calculations
switch funcformat
    
    % return cumulative distribution of wave troughs conditional on Hm0
    case 1
        H   = reshape(data.X  ,[],1);
        Hm0 = reshape(data.Hm0,1,[]);
        
		k   = 2*pi./lambda;
        nb1 = cosh(k.*d.*(2+cosh(2.*k.*d)));
        nb2 = 2.*(sinh(k.*d)).^3;
        nb3 = -1./(sinh(2.*k.*d));
		Rstar = k.*Hm0.*(nb1./nb2+nb3); %wave effective steepness
        Rstar = repmat(Rstar,length(H),1); %wave effective steepness
        output.CDF = 1 - exp(-8.*((Rstar.^2).^-1).*(sqrt(1-2.*Rstar.*(H*(Hm0.^-1)))-1).^2);
    
    % return cumulative distribution of Hmax conditional on Tmp % correct
    % or not??
    case 3
        Hmx  = reshape(data.Xmx,[],1);
        H    = reshape(data.X  ,1,[]);
        N    = data.N;
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^2-1)));
        
    otherwise
        disp('Error in call to Rayleigh T distribution function')
        return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Forristall_H(data,funcformat)
%
% Syntax: CDF = Forristall_H(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave heights 
%                 conditional on Hm0
%
% funcformat = 2: Return sample of maximum wave heigth in sea state
%
% funcformat = 3: Return cumulative distribution of Hmax
%                 conditional on Hmp
%
%
%

% default parameters
%
alfa = 2.126;
beta = 8.42; % note that this parameter is sometimes given as (beta^(1/alfa)/4) = 0.681
%
% perform calculations
switch funcformat
    
    % return cumulative distribution of wave heights conditional on Hm0
    case 1
        H          = reshape(data.X  ,[],1);
        Hm0        = reshape(data.Hm0,1,[]);
        output.CDF = 1 - exp(-1/beta*(4*H*(Hm0.^-1)).^alfa);

    % return random sample of max wave in sea state
    case 2
        if ~isfield(data,'P')
            data.P = rand(size(data.Hm0));
        end
        output     = (beta^(1/alfa)/4).*data.Hm0.*(-log(1-data.P.^(data.N.^-1))).^(1./alfa);
        
    % return cumulative distribution of Hmax conditional on Hmp
    case 3
        Hmx    = reshape(data.Xmx,[],1);
        H      = reshape(data.X  ,1,[]);
        N      = data.N;
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^alfa-1)));

    otherwise
        disp('Error in call to Forristall distribution function')
        return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Naess_H(data,funcformat)
%
% Syntax: CDF = Naess_H(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave heights 
%                 conditional on Hm0
%
% funcformat = 2: Return sample of maximum wave heigth in sea state
%
% funcformat = 3: Return cumulative distribution of Hmax
%                 conditional on Hmp
%
%
%
% perform calculations
switch funcformat
    
    % return cumulative distribution of wave heights conditional on Hm0
    case 1
        
        rho        = data.rho;
        alfaH      = 0.5*sqrt(1-rho);
        H          = reshape(data.X  ,[],1);
        Hm0        = reshape(data.Hm0,1,[]);
        output.CDF = 1 - exp(-(H*((alfaH*Hm0).^-1)).^2);
    
    % return random sample of max wave in sea state
    case 2
        rho        = data.rho;
        alfaH      = 0.5*sqrt(1-rho);
        if ~isfield(data,'P')
            data.P = rand(size(data.Hm0));
        end
        output     = alfaH.*data.Hm0.*(-log(1-data.P.^(data.N.^-1))).^(1/2);
        
    % return cumulative distribution of Hmax conditional on Hmp
    case 3
        
        Hmx    = reshape(data.Xmx,[],1);
        H      = reshape(data.X  ,1,[]);
        N      = data.N;
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^2-1)));

    otherwise
        disp('Error in call to N?ss distribution function')
        return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Battjes_Groenendijk_H(data,funcformat)
%
% 
HTRT = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 ...
    1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 ...
    2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.65 2.7 2.75 2.8 2.85 2.9 2.95 3];

H1T = [12.193 7.003 5.063 4.022 3.365 2.908 2.571 2.311 2.104 1.936 1.796 1.678 1.578 1.492 ...
    1.419 1.356 1.302 1.256 1.216 1.182 1.153 1.128 1.108 1.09 1.075 1.063 1.052 1.043 1.036 ...
    1.03 1.024 1.02 1.016 1.013 1.011 1.009 1.007 1.006 1.004 1.004 1.003 1.002 1.002 1.001 ...
    1.001 1.001 1.001 1 1 1 1 1 1 1 1 1 1 1 1 1];

H2T = [1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.06 1.061 1.061 1.062 1.064 1.066 1.069 1.073 ...
    1.077 1.083 1.09 1.097 1.106 1.116 1.126 1.138 1.15 1.162 1.175 1.189 1.203 1.217 1.231 ...
    1.246 1.261 1.275 1.29 1.305 1.32 1.334 1.349 1.363 1.378 1.392 1.407 1.421 1.435 1.449 ...
    1.462 1.476 1.49 1.503 1.516 1.529 1.542 1.555 1.568 1.58 1.593 1.605 1.617 1.63];
%
d     = reshape(data.d  ,1,[]);
alpha = data.slope;
%
switch funcformat

    case 1
        H    = reshape(data.X  ,[],1);
        Hm0  = reshape(data.Hm0,1,[]);
        res  = NaN(length(H),length(Hm0));
        
        Hm0mx = 0;
        for i=1:length(Hm0)
            Htr   = (0.35+5.8*tan(alpha))*d(i);
            Hrms  = (2.69+3.24*Hm0(i)/4/d(i))*Hm0(i)/4;
            Htrt  = Htr/Hrms;
            if Htrt < 0.05
                error('Hm0 unrealistically high for this water depth')
            end
            if Htrt > 3 % P(H|Hm0) follows Rayleigh distribution
                H1  = Hrms; 
                H2  = 1;
                Htr = max(H);
            else
                H1 = interp1(HTRT,H1T,Htrt)*Hrms;
                H2 = interp1(HTRT,H2T,Htrt)*Hrms;
            end
            % Return Htr at peak Hm0 
            if Hm0(i)>Hm0mx
                Hm0mx = Hm0(i);
                output.Htr = Htr;
            end

            res(:,i) = [1-exp(-(H(H<=Htr)/H1).^2) ; 1-exp(-(H(H>Htr)/H2).^3.6)];
        end
        output.CDF = res;

    case 2
        error('funcformat = 2 needs implementing for B&G distr.')
        
    case 3
        Hmx    = reshape(data.Xmx,[],1);
        H      = reshape(data.X  ,1,[]);
        N      = data.N;        
        Htr    = (0.35+5.8*tan(alpha))*d;
        output = [];
        if sum(H<Htr)
            output = [output,exp(-exp(-log(N)*((Hmx*(H(H<Htr).^-1)).^2-1)))];
        end
        if sum(H>=Htr)
            output = [output,exp(-exp(-log(N)*((Hmx*(H(H>=Htr).^-1)).^3.6-1)))];
        end
        
    otherwise
        disp('Error in call to Battjes and Groenendijk distribution function')
        return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Glukhovskiy_H(data,funcformat)
%
% Syntax: CDF = Glukhovskiy_H(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave heights 
%                 conditional on Hm0
%
% funcformat = 2: Return sample of maximum wave heigth in sea state
%
% funcformat = 3: Return cumulative distribution of Hmax
%                 conditional on Hmp
%
% Glukhovskiy wave height distribution as outlined in Battjes and
% Groenendijk (2000) paper (see also Klopman 1996)
%
%
% perform calculations
d     = reshape(data.d  ,1,[]);
%
switch funcformat

    
    
    % return cumulative distribution of wave heights conditional on Hm0
    case 1

        Hrms       = sqrt(1/2)*reshape(data.Hm0,1,[]);
        kappa      = 2./(1-0.7*Hrms./d);
        A          = (gamma(2./kappa+1)).^(kappa/2);
        H          = reshape(data.X  ,[],1);
        output.CDF = 1 - exp(-repmat(A,length(H),1).*(H*(Hrms.^-1)).^repmat(kappa,length(H),1));
    
        % Return kappa at peak Hm0 
        [~,imx] = max(Hrms);
        output.kappa_peak = kappa(imx(1));

    % return random sample of max wave in sea state
    case 2
        
        Hrms       = sqrt(1/2)*data.Hm0;
        kappa      = 2./(1-0.7*Hrms./data.d);
        A          = (gamma(2./kappa+1)).^(kappa/2);
        if ~isfield(data,'P')
            data.P = rand(size(data.Hm0));
        end
        output     = Hrms.*(1./A.*(-log(1-data.P.^(data.N.^-1)))).^(1./kappa);
        
    % return cumulative distribution of Hmax conditional on Hmp
    case 3
        
        Hrms   = sqrt(1/2)*max(data.Hm0);
        Hmx    = reshape(data.Xmx,[],1);
        H      = reshape(data.X  ,1,[]);
        N      = data.N;
        kappa  = 2/(1-0.7*Hrms(1)/d);
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^kappa-1)));

    otherwise
        disp('Error in call to Glukhovskiy distribution function')
        return
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Weibull_H(data,funcformat)
%
% Syntax: CDF = Weibull_H(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave heights 
%                 conditional on Hm0
%
% funcformat = 2: Return sample of maximum wave heigth in sea state
%
% funcformat = 3: Return cumulative distribution of Hmax
%                 conditional on Hmp
%
%
%
% default parameters
%
alfa = data.alpha; % this is the definition of alpha and beta used in Forristall 1978
beta = data.beta;

%
% perform calculations
switch funcformat
    
    % return cumulative distribution of wave heights conditional on Hm0
    case 1
        
        H          = reshape(data.X  ,[],1);
        Hm0        = reshape(data.Hm0,1,[]);
        output.CDF = 1 - exp(-1/beta*(4*H*(Hm0.^-1)).^alfa);

    % return random sample of max wave in sea state
    case 2
        if ~isfield(data,'P')
            data.P = rand(size(data.Hm0));
        end
        output     = (beta^(1/alfa)/4).*data.Hm0.*(-log(1-data.P.^(data.N.^-1))).^(1./alfa);

    % return cumulative distribution of Hmax conditional on Hmp
    case 3
        
        Hmx    = reshape(data.Xmx,[],1);
        H      = reshape(data.X  ,1,[]);
        N      = data.N;
        output = exp(-exp(-log(N)*((Hmx*(H.^-1)).^alfa-1)));

    otherwise
        disp('Error in call to Weibull_H distribution function')
        return
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[output] = Forristall_C(data,funcformat)
%
% Syntax: CDF = Forristall_C(data,funcformat)
%
% funcformat = 1: Return cumulative distribution of wave heights 
%                 conditional on Hm0
%
% funcformat = 2: Return sample of maximum crest heigth in sea state
%
% funcformat = 3: Return cumulative distribution of Hmax
%                 conditional on Hmp
%
%
%

% initialize
alfa = []; % Note that the definition of alpha and beta in Forristall 2000
beta = []; % is different from the definition used in Forristall 1978

% perform calculations
switch funcformat
    
    % return cumulative distribution of crests conditional on Hm0
    case 1

        C    = reshape(data.X  ,[],1);
        Hm0  = reshape(data.Hm0,1,[]);
        T01  = reshape(data.T01,1,[]);
        d    = reshape(data.d  ,1,[]);
        dim  = data.dim;
        alfa = data.alfa;
        beta = data.beta;
        res  = NaN(length(C),length(Hm0));
        
        if isempty(alfa) && isempty(beta)
            calcparm=1;
        else
            calcparm=0;
        end
        
        Hm0mx = 0;
        for i=1:length(Hm0)
        
            switch calcparm
                
                case 0  % alfa and beta have been specified in input. Do not recalc alfa and beta
                    output.alfa_peak = alfa;
                    output.beta_peak = beta;

                case 1 % Calculate alfa and beta

%                     L1 = DispRel(d(i),'T',T01(i));
                    L1 = m_DispRel(d(i),'T',T01(i)); % PDG: replaces the above line
                    k1 = 2*pi/L1;
                    S1 = 2*pi/9.81*Hm0(i)/T01(i)^2;
                    Ur = Hm0(i)/k1^2/d(i)^3;

                    if (isnumeric(dim) && dim==2) || strcmpi(dim,'2d')
                        alfa = 0.3536 + 0.2892*S1 + 0.1060*Ur;
                        beta = 2 - 2.1597*S1 + 0.0968*Ur^2;
                    elseif (isnumeric(dim) && dim==3) || strcmpi(dim,'3d')
                        alfa = 0.3536 + 0.2568*S1 + 0.0800*Ur;
                        beta = 2 - 1.7912*S1 - 0.5302*Ur + 0.2824*Ur^2; % PDG: 0.2824 was earlier mistakenly set to 0.284
                    end

                    % Return alfa and beta at peak Hm0 
                    if Hm0(i)>Hm0mx
                        Hm0mx = Hm0(i);
                        output.alfa_peak = alfa;
                        output.beta_peak = beta;
                    end

            end
            res(:,i) = 1 - exp(-(C/alfa/Hm0(i)).^beta);
        end
        output.CDF = res;

    % return random sample of max wave in sea state
    case 2
        
        if ~isfield(data,{'alfa','beta'})
            if ~isfield(data,'k1')
                data.k1 = w2k(2*pi./data.T01,[],data.d,9.81);
            end
            S1 = 2*pi/9.81*data.Hm0./data.T01.^2;
            Ur = data.Hm0./data.k1.^2./data.d.^3;
            if ~isfield(data,'dim') || data.dim==3 %3D
                alfa = 0.3536 + 0.2568*S1 + 0.0800*Ur;
                beta = 2 - 1.7912*S1 - 0.5302*Ur + 0.2824*Ur.^2; % PDG: 0.2824 was earlier mistakenly set to 0.284
            elseif data.dim==2 % 2D
                alfa = 0.3536 + 0.2892*S1 + 0.1060*Ur;
                beta = 2 - 2.1597*S1 + 0.0968*Ur.^2;
            end
        end
        if ~isfield(data,'P')
            data.P = rand(size(data.Hm0));
        end
        output = alfa.*data.Hm0.*(-log(1-data.P.^(data.N.^-1))).^(1./beta);        
        
    % return cumulative distribution of Cmax conditional on Cmp
    case 3
        
        Cmx    = reshape(data.Xmx,[],1);
        C      = reshape(data.X  ,1,[]);
        N      = reshape(data.N  ,1,[]);
        beta   = data.beta;
        output = exp(-exp(-log(N)*((Cmx*(C.^-1)).^beta-1)));
        
    otherwise
        disp('Error in call to Forristall distribution function')
        return
end

end
