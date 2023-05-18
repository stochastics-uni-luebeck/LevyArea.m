function LevyArea_setup
%LEVYAREA_SETUP Set up search paths for LevyArea.m
%   Run this function to setup the LevyArea.m package.
%
% See also: ADDPATH, RMPATH, SAVEPATH.

addpath(genpath(fileparts(mfilename('fullpath'))));

if savepath
    warning("Savepath failed. You may have to re-run LevyArea_setup" + ...
        " in each new MATLAB session.")
else
    disp("LevyArea.m was successfully added to the MATLAB path.")
end

end
