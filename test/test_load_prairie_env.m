
%%
% env_path = fullfile(env_dir, "Tseries_VivekNuria_holo_all.env")
% env_path = fullfile(env_dir, "Tseries_VivekNuria_holo_4.env")
% env_path = fullfile(env_dir, "Tseries_VivekNuria_15.env")
% env_path = fullfile(env_dir, "Tseries_VivekNuria_40.env")
%% prepare SEQUENTIAL HOLO STIM of individual neurons
%Load environment to prairie
% load environment
env_path = fullfile(env_dir, "Tseries_VivekNuria_holo_all.env")
exist(env_path)
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%Reps: 3200


%%
%Load environment to prairie
% load environment
env_path = fullfile(env_dir, "Tseries_VivekNuria_holo_4.env")
exist(env_path)
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()
%Reps: 1500
%%
%Load environment to prairie
% load environment
env_path = fullfile(env_dir, "Tseries_VivekNuria_40.env")
exist(env_path)
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()
%Reps: 75600
%%
%Load environment to prairie
% load environment
env_path = fullfile(env_dir, "Tseries_VivekNuria_15.env")
exist(env_path)
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()
%Reps: 27000
