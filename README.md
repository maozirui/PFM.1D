# PFM.1D
Phase Field modeling in 1D case

!!!!! IMPORTANT
 <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

By using this freeware, you are agree to the following:    
   1. you are free to copy and redistribute the material in any medium or format.
   2. you are free to remix, transform, and build upon the material for any purpose, even commercially.
   3. you must provide the name of the creator and attribution parties, a copyright notice, a license notice, a  
      disclaimer notice, and a link to the material (https://github.com/maozirui/PFM.1D).
   4. users are entirely at their own risk using this freeware and techniques. 
 
 Before use, please read the License carefully.

 
Introduction:   
This script is developed to simulate the phase field evolution phenomenon governed by the Allen-Cahn equation and Cahn-Hillard equation in 1D case. The attached is Matlab GUI freeware. 

Author Information:   
Zirui Mao | Department of Material Science and Engineering | TAMU
Date last modified: March., 10th, 2020

How to use this freeware:   
1. Download the documents 'Diffusion.m' and 'Diffusion.fig'; 
2. Open the 'Difussion.m' with Matlab.              
3. Before running it, please ensure the two documents locate in the same folder.                                               
4. Click 'Run' in Matlab toolstrip under 'EDITOR'. 
5. The GUI (User Interface) shall show up.
6. Specify the values for each user-define parameters;     
7. Click the 'Run' botton in GUI when everything is ready;  
8. The instant result will be plotted in the right panel;   
9. Click the 'Stop' button in GUI to stop the calculation at any time;
10. An animation will be saved in the document-located folder after each time of calculation. Please rename it immediately otherwiseit will be overwritten.


% %%%%%%%%% Numerical model:                                             
%  /|\                             |                                     
%   |                              |                                     
%   |                              |                                     
%   |                                                                    
%   |                   ___________C2___________                         
%   |                  |           |            |                        
%   |                  |           |            |                        
%   |                  |                        |                        
%   |_______C1_________|<-------- w=r*L ------->|__________C1________    
%   |                                                                    
%   |<---------------------------- L ------------------------------->|   %
%   |__________________________________________________________________\ %
%   0                              |                                   / %
%                                  |                                     %
% %%%%%%%%% Governing equation:                                          %
%                                                                        %
%   Allen-Cahn equation:                                                 %
%   dc/dt = - M * [\Deltaf * (df/dc) - \gamma * div(c)]                  %
%   Cahn-Hillar equation:                                                %
%   dc/dt = M * div(\Deltaf * (df/dc) - \gamma * div(c))                 %
%                                                                        %
% where c = phase field (also always indicated by \phi) c in [-1, 1]     %
%       M = Mobility coefficient                                         %
%       f = local potential energy = (c^2-1)^2                           %
%       \Deltaf = local potential energy difference                      %
%       \gamma = gradient energy coefficient                             %
%       L = length of 1D domain                                          %
%       w = width of the central segment controlled by the factor f      %
%       r = ratio of the central segment to L.  Note: r belong (0,1)     %
%      C1 = one constant defining the initial c in both sides            %
%      C2 = another constant defining the initial c in central segment   %
%       iter = total number of iterations                                %
%                                                                        %
% %%%%%%%%% Boundary Conditions:                                         %
% Periodical boundary is applied to the both ending nodes, i.e.,         %
%      x(N) | x(1) x(2) .................. x(N-1) x(N) | x(1)            %
%                                                                        %
% %%%%%%%%% Finite Difference Approximation:                             %
% The second-order derivative of c in the govering equation (1) is       %
% approximated with the 2nd order accurate central Finite Difference     %
% scheme, i.e.,                                                          %
%                                                                        %
%  d^2(c)   c(i+1) - 2c(i) + c(i-1)                                      %
%  ------ = -----------------------                                      %
%   dx^2           (dx)^2                                                %
%                                                                        %
% %%%%%%%%% inputs:                                                      %
% L, r, C1, C2; M, \gamma, \Deltaf, iter;  A-C or C-H                    %
%                                                                        %  
