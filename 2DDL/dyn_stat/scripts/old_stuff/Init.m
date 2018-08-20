classdef Init
    %INIT initializes the script according to standard parameters such as
    %clearing the variables, clearing the figures and the command window
    %   2 methods: 1) clear(), which holds clear vars close all and clc; 2)
    %   print_setup() which sets print properties suggested by an IEEE
    %   paper. 
    
    properties
    end
    
    methods (Static)
        function clear()
            clear all; 
            close all; 
            clc;
        end
        function print_setup()
            set(0,'DefaultTextFontName','Times',...
                'DefaultTextFontSize',18,...
                'DefaultAxesFontName','Times',...
                'DefaultAxesFontSize',18,...
                'DefaultLineLineWidth',1,...
                'DefaultLineMarkerSize',7.75,...
                'defaultfigurecolor','w');
            set(groot, 'defaultAxesTickLabelInterpreter','latex');...
            %set(groot, 'defaultLegendInterpreter','latex');
            %set(groot, 'DefaultTextInterpreter','latex');
%             set(0,'defaultfigurecolor',[.4 .1 .5])
        end
    end
    
end

