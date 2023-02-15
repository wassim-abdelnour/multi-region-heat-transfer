#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Company : Temisth
# Creation Date: 25 - Janvier - 2023
# version ='2.0'
# ---------------------------------------------------------------------------
""" Cold Plates Thermo-Physical Caracterization using openFoam """
# ---------------------------------------------------------------------------
#import math
#import numpy as np
#import matplotlib
#import cmath
#import gmsh
#import os
# ---------------------------------------------------------------------------
import foamControl
import foamSc
import foamGenMesh
import foamProcess
# ---------------------------------------------------------------------------


'-------------------------- Temisth --------------------------'
'***Technology Leader for High Performances Thermal systems***'
'-------------------------------------------------------------'

' Cold Plates Thermo-Physical Caracterization using openFoam '

'-------------------------------------------------------------'

'version 4.0 - 6 Fev 2023'

'-------------------------------------------------------------'

#simulation_name, plate_width, plate_height, plate_thickness, longueur, largeur, spacing, nombre_de_passes, perc_of_r, debit_massique, rho, mu, Cp_f, kf, Q, Tin, rho_s, Cp_s, ks, number_of_iterations, Px, Py, Pz = foamControl.get_inputs(sim_log)

plate_width = 0.095 #0.095
plate_height = 0.2 #0.2
plate_thickness = 0.014 #0.014
longueur = 0.01 #0.01
largeur = 0.01 #O.01
spacing = 0.02 #0.02
nombre_de_passes = 1 #1
perc_of_r = 100 #100
rho = 1044 #1044
mu = 0.00157 #0.00157
Cp_f = 3680 #3680
kf = 0.473 #0.473
Q = 150 #150
Tin = 273 #273
rho_s = 2698.9 #2698.9
Cp_s = 900 #900
ks = 210 #210
pack_of_iterations = 100 #multiple de 10
convergence_criterion = 1 #5%
Px = 3 #8
Py = 3 #8
Pz = 1 #1

debit_volumique = [x*0.5 for x in range(2, 21, 1)] #l/s
debit_massique_iter = [x*rho / 1000 / 60 for x in (debit_volumique)] #kg/s

i = 2
for debit_massique in debit_massique_iter:
    #simulation header
    simulation_name = foamControl.get_simulation_name(i, debit_massique)
    sim_log = open(simulation_name + "/gen.log", "w")
    foamControl.print_simulation_recap(plate_width, plate_height, plate_thickness, longueur, largeur, spacing, perc_of_r, nombre_de_passes, debit_massique, rho, mu, Tin, Cp_f, kf, Q, rho_s, Cp_s, ks, Px, Py, Pz, sim_log)

    #compute turbulence parameters
    k, eps, nut, omega, l = foamControl.compute_turbulence_parameters(longueur, largeur, debit_massique, rho, mu, sim_log)

    #create simulation setup
    foamSc.create_openFoam_simulation_setup(simulation_name, debit_massique, rho, mu, Cp_f, kf, Q, Tin, rho_s, Cp_s, ks, Px, Py, Pz, k, eps, nut, omega, l, pack_of_iterations, pack_of_iterations)

    #generate meshes
    foamGenMesh.generate_fluid_mesh(plate_width, plate_height, plate_thickness, longueur, largeur, spacing, nombre_de_passes, perc_of_r, debit_massique, rho, mu, sim_log)
    foamGenMesh.generate_solid_mesh(plate_width, plate_height, plate_thickness, longueur, largeur, spacing, nombre_de_passes, perc_of_r, sim_log)
    foamGenMesh.Convert_Mesh_to_Foam_Format(simulation_name, sim_log)

    #lauch the simulation
    number_of_iterations = foamControl.run_simulation(simulation_name, pack_of_iterations, Q, debit_massique, Cp_f, convergence_criterion, Px, Py, Pz, sim_log)

    #post process
    foamProcess.compute_head_losses(simulation_name, number_of_iterations, Px, Py, Pz, sim_log)
    foamProcess.compute_thermal_resistance(simulation_name, number_of_iterations, Q, Px, Py, Pz, sim_log)
    foamProcess.compute_Yplus(simulation_name, number_of_iterations, Px, Py, Pz, sim_log)
    foamProcess.plot_convergence_curve(simulation_name, number_of_iterations, pack_of_iterations, Q, rho, Cp_f, debit_massique, Tin, sim_log)

    #for visualization
    foamControl.reconstruct_simulation(simulation_name, sim_log)
    foamControl.Extract_Solver_Details(simulation_name, number_of_iterations, pack_of_iterations, sim_log)
    foamProcess.plot_residuals(simulation_name, number_of_iterations, pack_of_iterations, sim_log)
    foamProcess.compute_simulation_time(simulation_name, number_of_iterations, pack_of_iterations, sim_log)

    #end of prog
    sim_log.close()
    i = i+1
print()
#end_of_prog = input('Press Enter to exit')
