%% Script to calculate the gradiant of the stat_wS limits 
clear all;close all;clc;


syms L l ay cy 

A2 = ay*(cy-L)+2*L*cy-l*(ay+L);
A3 = -ay*(cy-L)-2*L*cy-l*(ay+L);
B2 = ay-l;
B3 = -(ay+l);

eq2 = L-(A2/B2);
eq3 = L-(A3/B3);

gradeq2 = gradient(eq2,[ay,cy]);

gradeq3 = gradient(eq3,[ay,cy]);
