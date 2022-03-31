    clear
    clc
% %=======================================================Background paths
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/core')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/solve')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/exp')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/cross')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/fmex')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/misc')
%     addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/solve')
%=======================================================================
%---Load TT Toolbox directories and its subfolder (copy content from
%TT-Tool
  run ('/Users/trangdinh/Documents/Cod.fsp-qtt-sw/TT-Toolbox/setup.m')
%---Load characteristics for the specific model
    module_model_Michaelis_Menten
    % module_model_Michaelis_Menten_big
    % module_model_Gene_toggle
    % module_model_p53
%---Load parameters for the solvers
    module_global
%---Run the program
    example
