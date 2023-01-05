%"""
%###############################################################################,
%# Created Date:    '2022-06-21'                                                #,
%# Author(s):       Mohammad Madani - mohm@dhigroup.com                         #,
%#                  Benjamin Hernandez Alfaro - beha@dhigroup.com               #,
%# Encoding:        UTF-8                                                       #,
%# Language:        MATLAB                                                      #,
%# ---------------------------------------------------------------------------- #,
%# This script is a template to perform EVA                                     #,
%# sensitivity = 1 if you want to perform sensitivity                           #,
%#                                                                              #,
%# ---------------------------------------------------------------------------- #,
%# Copyright (c) (2022) DHI Water & Environment, Inc.                           #,
%################################################################################
%"""


load cmp_data.mat


M_CS = m_structure(name , xyz , ttt_model , legend_model , fname_model ,'Hm0' , colums_model(1), bins_hm0);
O_CS = m_structure(name , xyz , ttt_obs , legend_obs , fname_obs ,'Hm0' , colums_obs(1), bins_hm0);


M_CS.fontsize = 6;

m_timeseries(M_CS)
m_compare(M_CS, O_CS,'tmin',datenum([2019 12 01 00 00 00]),'tmax',datenum([2020 07 01 00 00 00]))
m_scatter(M_CS, O_CS)



