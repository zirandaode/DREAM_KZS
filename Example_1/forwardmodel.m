function [SimRR] = forwardmodel(x)
% Runs the HYMOD model

% Define local variables
persistent MaxT PET Precip F

% Load the boundary conditions only once
if isempty(F)
    
    % Load the Leaf River data        
    load bound.txt;
    % Only use two years of data
    MaxT = 795;
    % Extract the PET
    PET = bound(1:MaxT,5); 
    % Extract the precipitation
    Precip = sum(bound(1:MaxT,6:9),2);
    % Area factor to translate HYMOD output in mm/d to m3/s (calibration data); (area Leaf River is 1944 km2)
    F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);

end;

% -------------------------------------------------------------------------
%                           Model script
% -------------------------------------------------------------------------

% Define the parameters
cmax = x(1); bexp = x(2); alpha = x(3); Rs = x(4); Rq = x(5);

% Set the initial states
x_loss = 0.0;

% Initialize slow tank state
x_slow = 0; % --> works ok if calibration data starts with low discharge

% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0; outflow = [];

% Now loop over the forcing data
for t = 1:MaxT
    
    % Assign precipitation and evapotranspiration
    Pval = Precip(t,1); PETval = PET(t,1);
    
    % Compute excess precipitation and evaporation
    [ER1,ER2,x_loss] = excess(x_loss,cmax,bexp,Pval,PETval);
    
    % Calculate total effective rainfall
    ET = ER1 + ER2;
    
    % Now partition ER between quick and slow flow reservoirs
    UQ = alpha*ET; US = (1-alpha)*ET;
  
    % Route slow flow component with single linear reservoir
    [x_slow,QS] = linres(x_slow,US,Rs);
    
    % Route quick flow component with linear reservoirs
    inflow = UQ; 
    
    for k = 1:3
        % Linear reservoir
        [x_quick(k),outflow] = linres(x_quick(k),inflow,Rq); inflow = outflow;
    end;

    % Compute total flow for timestep
    output(t,1) = (QS + outflow); %#ok<*AGROW>
    
end;

SimRR = F * output(65:MaxT,1);

% -------------------------------------------------------------------------
function [ER1,ER2,xn] = excess(x_loss,cmax,bexp,Pval,PETval)
% this function calculates excess precipitation and evaporation

xn_prev = x_loss;
ct_prev = cmax*(1-power((1-((bexp+1)*(xn_prev)/cmax)),(1/(bexp+1))));
% Calculate Effective rainfall 1
ER1 = max((Pval-cmax+ct_prev),0.0);
Pval = Pval-ER1;
dummy = min(((ct_prev+Pval)/cmax),1);
xn = (cmax/(bexp+1))*(1-power((1-dummy),(bexp+1)));
% Calculate Effective rainfall 2
ER2 = max(Pval-(xn-xn_prev),0);

% Alternative approach
evap = (1-(((cmax/(bexp+1))-xn)/(cmax/(bexp+1))))*PETval; % actual ET is linearly related to the soil moisture state
xn = max(xn-evap, 0); % update state

%evap = min(xn,PETval);
%xn = xn-evap;

% -------------------------------------------------------------------------

function [x_slow,outflow] = linres(x_slow,inflow,Rs)
% Linear reservoir
x_slow = (1-Rs)*x_slow + (1-Rs)*inflow;
outflow = (Rs/(1-Rs))*x_slow;
