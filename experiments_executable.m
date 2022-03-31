clear
clc
%=======================================================Background paths
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/core')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/solve')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/exp')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/cross')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/fmex')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/misc')
addpath('/Users/TrangDinh/Google Drive/Dinh Tensor Project/Tensor Toolbox MATLAB/TT-Toolbox-Trang/solve')
%=======================================================================
% %---Load mat file for the SSA when run the specific model the 1st time
%     Note: Comment the below if run the model the 1st time
%     and Uncomment them the below if..end after the 1st run

% load data_SSA_goutsias.mat
% load data_SSA_goutsias_10000_200603_2.mat
% load data_SSA_goutsias_S100000.mat
% load data_goutsias_SSA_10000_init_[3 3 3 3 3 3]_tf_100_200604_1.mat

load data_SSA_goutsias_10000_init_[443433]_tf_100_200604_1.mat
% %---Load mat file for the specific model
load data_goutsias.mat
% load data_goutsias_10000_200603_2.mat
% load data_goutsias_S10_HuyFalse.mat
%---Load file to output plots and internal/external statistics
plot_compact
