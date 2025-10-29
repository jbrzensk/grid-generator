% Plot grids for Aneesh to solve
%
%
close all; clear all; clc;

nptsx = 81;
nptsy = 41;

[ Xa, Ya ] = makeGridA( nptsx, nptsy );

[ Xc, Yc ] = makeGridC( nptsx, nptsy );

[ Xh, Yh ] = makeGridHardC( nptsx, nptsy );
