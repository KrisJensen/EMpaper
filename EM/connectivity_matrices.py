#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:07:30 2020

@author: kris
"""

from buildCX import buildCX, connectivity

show = True #show the neurons as we import them
vol = buildCX(show=show) #construct a 'volume' instance from the neurons

connectivity(vol, show) #get connectivity matrices

