[fileName, pathName] = uigetfile('*.csv');
data = csvread(fullfile(pathName,fileName));