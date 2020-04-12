clear all
close all
clc

load stress1.out
load strain1.out
load stress2.out
load strain2.out

plot(strain1(:,2),stress1(:,2),'b.')
hold on
%plot(strain2(:,2),stress2(:,2),'r.')