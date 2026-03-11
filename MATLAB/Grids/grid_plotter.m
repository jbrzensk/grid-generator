% Plot grids for Aneesh to solve
%
%
close all; clear all; clc;

nptsx = 32;
nptsy = 21;

[ Xa, Ya ] = makeGridA( nptsx, nptsy );

[ Xc, Yc ] = makeGridC( nptsx, nptsy );

[ Xh, Yh ] = makeGridHardC( nptsx, nptsy );

[ Xv, Yv ] = makeGridChevron( nptsx, nptsy );

[ Xf, Yf ] = makeGridFish( nptsx, nptsy );
