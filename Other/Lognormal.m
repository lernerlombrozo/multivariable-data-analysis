clear; clc;

x =[4; 8; 16; 30; 50];
y=[25; 125; 100; 75; 30; 5];
ycum=[25; 150; 250; 325; 355; 360]
mu=11;
sigma=2;
crit =lognpdf(x,mu,sigma)
