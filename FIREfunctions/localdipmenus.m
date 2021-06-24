% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Oct 2012

function [menulist,excludelist] = localdipmenus(menulist)
I = size(menulist,1)+1;
menulist{I,1} = 'FIRE';
menulist{I,2} = {'fire_ims','fire_locs','-','frc','flc','fpc','-'};
excludelist = {'isect','sphere_tesselation','wkk','gfca3D_sub','frctoresolution'};