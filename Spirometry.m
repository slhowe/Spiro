clc
close all
clear
debug_on_warning;
debug_on_error;

load SpirometryData.mat;
loops = data.Loops;
band = data.Banding;
normal = data.Normal;
inflated = data.Inflated;

% start normal
flow = normal.flow;
pressure = normal.pressure;
