clear; clc; close all;

param = 0;

disp('Creating Masks');
tic;
if ~exist('Files/l0.mat', 'file')
    createPatchMasks(param);
end

disp('Building Tree');

dirInfo = dir('Images/');

disp('PNG files in the directory:');
for i = 1:length(dirInfo)
    if ~dirInfo(i).isdir && endsWith(dirInfo(i).name, '.png', 'IgnoreCase', true)
        fprintf('%d) %s\n', i-2, dirInfo(i).name);
    end
end

menu = input("select image: ");
switch num2str(menu)
    case "1"
        I = im2double(imread('Images/CC.png'));
    case "2"
        I = im2double(imread('Images/Sines.png'));
    case "3"
        I = im2double(imread('Images/Sqr.png'));
end


fastRun = false;
if fastRun
    I = imresize(I,[65 65]);
end

sigma = 0.1;
[res,im] = runIm(I,false,param);
disp('Displaying Edges');
figure; 
subplot(1,2,1);imshow(res);
subplot(1,2,2);imshow(I);
toc;