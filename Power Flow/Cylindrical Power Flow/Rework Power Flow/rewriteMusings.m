%% write simpler code for power flow stuff...

%% break it up into more classes?
% Lets start high level and write test proceedure and then work down


%%%% Basic Functions
%% read in data from SLDV
a SLDV_parser obj takes in the SLDV ASCII File(s) and outputs an SLDV_data object.
This object has values of frequency coordinate displacement and velocity
%% compute surfaces from data points
A surface creator object takes in the SLDV_data object and using splines out puts a surface object
%% take derivatives of surfaces
The surface spline is fed into a derivative generator and it can output requested derivatives of the surface

%% compute resultants
A power flow resultant object uses the surfaces and derivative generators to compute power

%% visualize power flow
The pf_visualizer takes in power flow object and can then visualize the results. 