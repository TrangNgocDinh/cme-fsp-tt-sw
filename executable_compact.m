clear
clc
%=======================================================Background paths
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/core')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/solve')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/exp')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/cross')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/fmex')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/misc')
addpath('/Users/TrangDinh/GoogleDrive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/solve')
%=======================================================================
%---Load characteristics for the specific model
% module_model_Michaelis_Menten
% module_model_Michaelis_Menten_big
% module_model_Gene_toggle
module_model_p53
% module_model_Goutsias
%---Load parameters for the solvers
module_global
%---Run the program
example_compact
%---Run the experiments program
%   experiments_executable
