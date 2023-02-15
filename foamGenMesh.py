import math
import numpy as np
import gmsh
import os

def print_GenMesh_log(print_GenMesh_log, sim_log):
    if(print_GenMesh_log == "\n"):
        print()
        sim_log.write(print_GenMesh_log)
    else:
        print(print_GenMesh_log)
        sim_log.write(print_GenMesh_log + "\n")

list_of_walls = []
list_of_volumes = []
def update_list_of_walls_and_volumes(Surfaces):
    for i in range (1, 24, 6):
        list_of_volumes.append(Surfaces[i][1])
        for j in range(i+1, i+5):
            list_of_walls.append(Surfaces[j][1]) #index of external curve in curve_loops

list_of_solid_walls = []
Surface_Loop_Faces = []
def update_solid_walls(Surfaces):
    for i in range(1, 14, 4):
        list_of_solid_walls.append(Surfaces[i][1]) #index of external curve in curve_loops(2-3-4-5)

def generate_fluid_mesh(plate_width, plate_height, plate_thickness, longueur, largeur, spacing, nombre_de_passes, perc_of_r, debit_massique, rho, mu, sim_log):
#plate_width, plate_height, longueur, largeur, spacing, nombre_de_passes and perc_of_r are used to draw the duct
#plate thickness is not used, but we added it for consistency
#debit_massique, rho and mu are used to generate the mesh with convinient Y+
#sim_log is used to continue filling the simulation's general log file
    #plate_width = plate width
    #plate_height : plate height
    #plate_thickness : plate thickness
    #longueur : duct width
    #largeur : duct thickness
    #spacing : distance between duct's borders and plate's borders on each side (left, right, and top)
    #nombre_de_passes : number of turns the duct makes inside the plate
    #perc_of_r : steepness of the connexion's curvature (the lower is perc_of_r, the stepper is the curvature)
    #debit_massique : mass flow rate inside the duct
    #rho : density of the fluid
    #mu : dynamic viscosity of the fluid

    #inlet area = Pi*r2
    S = longueur * largeur
    P = 2 * (longueur + largeur)
    Hd = 4*S / P

    #debit_massique = rho * V * S
    Vmoyenne = debit_massique / (rho*S)

    #free-stream velocity = vitesse maximale d'un profile paraboloique au meme débit
    Vinf = 2*Vmoyenne #à revoir, puisque le tube n'est plus symétrique. Mais bonne approximation pour l'instent

    #The caracteristic length for Re calculation is the hydraulic diameter
    #Hydraulic Diameter = 4*Surface_Area/Perimeter
    #Hydraulic Number of a Rectangular cross section is 4*(L*l) / 2*(L+l)
    Re = rho * Vmoyenne * Hd / mu

    Yplus = 1 #constant since we use k-w SST
    element_size = longueur / 40 #will be automatically set based on inlet width

    block_width = longueur/2
    block_height = largeur/2

    #Cf = 0.078*Re^(0.25)
    Cf =0.078*(math.pow(Re, 0.25))

    #Yp = (mu/rho) * (1/sqrt(0.5*Cf*Vinf²)) * Yplus
    Yp = mu*Yplus*10 / (rho* math.sqrt(0.5*Cf*Vinf*Vinf))

    print_GenMesh_log("Target Y+ = " + str(Yplus), sim_log)
    print_GenMesh_log("Target Base Element Size = " + "{:e}".format(element_size), sim_log)

    print_GenMesh_log("Centre of neerest cell to the interface will be Yp = " + "{:e}".format(Yp) + " m\n", sim_log)

    first_element_of_expansion = 2*Yp
    ratio_lastElement_firstElement = element_size/(2*Yp)

    if (ratio_lastElement_firstElement > 1):
        geometric_exp_ratio_l = (block_width - first_element_of_expansion)/(block_width - element_size)
        geometric_exp_ratio_L = (block_height - first_element_of_expansion)/(block_height - element_size)
        number_of_cells_l = int( math.log(element_size/first_element_of_expansion)/math.log(geometric_exp_ratio_l))+1
        number_of_cells_L = int( math.log(element_size/first_element_of_expansion)/math.log(geometric_exp_ratio_L))+1

        print_GenMesh_log("number of cells in the horizontal direction : " + str(number_of_cells_l*2), sim_log)
        number_of_nodes_l = number_of_cells_l + 1
        print_GenMesh_log("number of nodes in the horizontal direction : " + str(number_of_nodes_l*2-1), sim_log)

        print_GenMesh_log("number of cells in the vertical direction : " + str(number_of_cells_L*2), sim_log)
        number_of_nodes_L = number_of_cells_L + 1
        print_GenMesh_log("number of nodes in the vertical direction : " + str(number_of_nodes_L*2-1), sim_log)

        coeff_l = np.zeros(number_of_cells_l)
        coeff_l += (2*Yp) #Cz Yp is the cell center and we are searching for the cell width
        coeff_l[number_of_cells_l-1] -= block_width
        coeff_l = np.roots(coeff_l)
        coeff_l = [x for x in coeff_l if (x.imag==0 and x>0 and x!=1)]
        #In Openfoam, the input is the ration between first and last element sizes (here the input is the geometric coefficient)    
        #    expansion_ratio = pow(coeff[0].real, number_of_cells-1)
        expansion_ratio_l = coeff_l[0].real
        print_GenMesh_log("\texpansion ratio of progression in the horizontal direcrion : " + "{:.2f}".format(expansion_ratio_l), sim_log)

        coeff_L = np.zeros(number_of_cells_L)
        coeff_L  += (2*Yp) #Cz Yp is the cell center and we are searching for the cell width
        coeff_L [number_of_cells_L-1] -= block_height
        coeff_L  = np.roots(coeff_L )
        coeff_L  = [x for x in coeff_L  if (x.imag==0 and x>0 and x!=1)]
        expansion_ratio_L = coeff_L [0].real
        print_GenMesh_log("\texpansion ratio of progression in the vertical direcrion : " + "{:.2f}".format(expansion_ratio_L), sim_log)
    else : 
        expansion_ratio_l = 1
        expansion_ratio_L = 1
        number_of_cells_l = max(int(block_width/(2*Yp)), 10)
        number_of_cells_L = max(int(block_height/(2*Yp)), int(block_height*10/block_width) )
        number_of_nodes_l = number_of_cells_l + 1
        number_of_nodes_L = number_of_cells_L + 1
        print_GenMesh_log("number of cells in the horizontal direction : " + str(number_of_cells_l*2), sim_log)
        print_GenMesh_log("number of cells in the vertical direction : " + str(number_of_cells_L*2), sim_log)
        print_GenMesh_log("no expansion nor progression", sim_log)

    ##############################################################

    space = plate_width - 2*spacing - 2*longueur
    height = plate_height - spacing - longueur

    R = perc_of_r * (longueur/2) /100

    xC1 = -space/2 #position de l'entrée
    h = height - longueur/2 - R #hauteur du canal d'entrée
    h_mid = h - (longueur/2+R) - longueur/2 -spacing #hauteur des canaux intermédiaires
    s = (space - 2*(longueur/2+R)*(2*nombre_de_passes-1)) / (2*nombre_de_passes-1) #espacements entre les canaux

    dz = int(h/element_size/6.5) +1 #nombre d'éléments dans les tubes d'entrée et de sortie
    dz_mid = int(h_mid/element_size/6.5) +1 #nombre d'éléments dans les tubes verticaux
    dzH = int(s/element_size/6.5) +1 #nombre d'éléments dans les tubes horizontaux
    dzC = int((math.pi*(longueur/2+R)) / (2*element_size*6.5)) +1 #nombre d'éléments dans les coins

    number_of_elements = number_of_cells_l * number_of_cells_L *(dz*2 + dzC*2 + (nombre_de_passes-1)*2*dz_mid + ((nombre_de_passes-1)*2+1)*dzH + ((nombre_de_passes-1)*4+2)*dzC)
    print_GenMesh_log("number of nodes in main ducts : " + str(dz), sim_log)
    print_GenMesh_log("number of nodes in auxiliary vertical ducts : " + str(dz_mid), sim_log)
    print_GenMesh_log("number of nodes in auxiliary horizontal ducts : " + str(dzH), sim_log)
    print_GenMesh_log("number of nodes in duct connectors : " + str(dzC), sim_log)
    print_GenMesh_log('\n', sim_log)
    print_GenMesh_log("total number of elements : " + str(number_of_elements), sim_log)
    print_GenMesh_log('\n', sim_log)

    print_GenMesh_log('drawing fluid geometry ...', sim_log)

    ###################Generation of fluid Mesh###################
    gmsh.initialize()

    #Centres des Cercles
    point1 = gmsh.model.geo.add_point(xC1, 0, 0)

    point2 = gmsh.model.geo.add_point(xC1-longueur/2, 0, -largeur/2)
    point3 = gmsh.model.geo.add_point(xC1+longueur/2, 0, -largeur/2)
    point4 = gmsh.model.geo.add_point(xC1+longueur/2, 0, +largeur/2)
    point5 = gmsh.model.geo.add_point(xC1-longueur/2, 0, +largeur/2)

    #Sommets du Volume Intérieur
    point6 = gmsh.model.geo.add_point(xC1-longueur/2, 0, 0)
    point7 = gmsh.model.geo.add_point(xC1+longueur/2, 0, 0)

    point8 = gmsh.model.geo.add_point(xC1, 0, -largeur/2)
    point9 = gmsh.model.geo.add_point(xC1, 0, +largeur/2)

    #Courbes du Volume Intérieur
    line = [0 for x in range(12)]
    line[0] = gmsh.model.geo.add_line(point2, point8)
    line[1] = gmsh.model.geo.add_line(point6, point1)
    line[2] = gmsh.model.geo.add_line(point5, point9)
    line[3] = gmsh.model.geo.add_line(point3, point8)
    line[4] = gmsh.model.geo.add_line(point7, point1)
    line[5] = gmsh.model.geo.add_line(point4, point9)

    line[6] = gmsh.model.geo.add_line(point5, point6)
    line[7] = gmsh.model.geo.add_line(point9, point1)
    line[8] = gmsh.model.geo.add_line(point4, point7)
    line[9] = gmsh.model.geo.add_line(point2, point6)
    line[10] = gmsh.model.geo.add_line(point8, point1)
    line[11] = gmsh.model.geo.add_line(point3, point7)


    # loops of base geometry
    curve_loop1 = gmsh.model.geo.add_curve_loop([ line[2], -line[6], -line[1], line[7] ])
    curve_loop2 = gmsh.model.geo.add_curve_loop([ line[5], -line[8], -line[4], line[7] ])
    curve_loop3 = gmsh.model.geo.add_curve_loop([ line[9], -line[0], -line[10], line[1] ])
    curve_loop4 = gmsh.model.geo.add_curve_loop([ line[3], -line[11], -line[4], line[10] ])

    # surfaces of base geometry
    face1 = gmsh.model.geo.add_plane_surface([curve_loop1])
    face2 = gmsh.model.geo.add_plane_surface([curve_loop2])
    face3 = gmsh.model.geo.add_plane_surface([curve_loop3])
    face4 = gmsh.model.geo.add_plane_surface([curve_loop4])

    for i in range(6):
        gmsh.model.geo.mesh.setTransfiniteCurve(line[i], number_of_nodes_l, "Progression", expansion_ratio_l)

    for i in range (6,12):
        gmsh.model.geo.mesh.setTransfiniteCurve(line[i], number_of_nodes_L, "Progression", expansion_ratio_L)

    gmsh.model.geo.mesh.setTransfiniteSurface(face1)
    gmsh.model.geo.mesh.setTransfiniteSurface(face2)
    gmsh.model.geo.mesh.setTransfiniteSurface(face3)
    gmsh.model.geo.mesh.setTransfiniteSurface(face4)

    gmsh.model.geo.mesh.setRecombine(2, face1)
    gmsh.model.geo.mesh.setRecombine(2, face2)
    gmsh.model.geo.mesh.setRecombine(2, face3)
    gmsh.model.geo.mesh.setRecombine(2, face4)

    #Surface1 = gmsh.model.geo.extrude([(dimension (always2, face_tag)], extrusion_x, extrusion_y, extrusion_z, [list of number of sub-divisions in each layer], [list of fraction of the layer in the block (empty if we have only one layer)], c)

    #Cylindre Vertical Gauche
    Surfaces = gmsh.model.geo.extrude([(2, face1), (2, face2), (2, face3), (2, face4)], 0, h, 0, [dz], [], True)
    update_list_of_walls_and_volumes(Surfaces)

    #print_GenMesh_log(Surfaces1[0])
    #attention : Surfaces1[0] est un vecteur, donc Surfaces1 est une matrice
    #Surfaces1 = gmsh.model.geo.revolve([(dimension (always2, face_tag)], rot_center_x, rot_center_y, rot_center_z, rot_axis_x, rot_axis_y, rot_axis_z, rot_angle, [list of number of sub-divisions in each layer], [list of fraction of the layer in the block (empty if we have only one layer)], True/False to recombine mesh)

    #Connexion Gauche
    Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], xC1+longueur/2+R, h, 0, 0, 0, -1, math.pi/2, [dzC], [], True)
    update_list_of_walls_and_volumes(Surfaces)

    #Cylindre Horizontal
    Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], s, 0, 0, [dzH], [], True)
    update_list_of_walls_and_volumes(Surfaces)

    for i in range(nombre_de_passes - 1):
        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], xC1+(1+4*i)*(longueur/2+R)+(1+2*i)*s, h, 0, 0, 0, -1, math.pi/2, [dzC], [], True)
        update_list_of_walls_and_volumes(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], 0, -h_mid, 0, [dz_mid], [], True)
        update_list_of_walls_and_volumes(Surfaces)

        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], xC1+(3+4*i)*(longueur/2+R)+(1+2*i)*s, height-h_mid-(longueur/2+R), 0, 0, 0, -1, -math.pi/2, [dzC], [], True)
        update_list_of_walls_and_volumes(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], s, 0, 0, [dzH], [], True)
        update_list_of_walls_and_volumes(Surfaces)

        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], xC1+(3+4*i)*(longueur/2+R)+(2+2*i)*s, height-h_mid-(longueur/2+R), 0, 0, 0, -1, -math.pi/2, [dzC], [], True)
        update_list_of_walls_and_volumes(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], 0, h_mid, 0, [dz_mid], [], True)
        update_list_of_walls_and_volumes(Surfaces)

        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], xC1+(5+4*i)*(longueur/2+R)+(2+2*i)*s, h, 0, 0, 0, -1, math.pi/2, [dzC], [], True)
        update_list_of_walls_and_volumes(Surfaces)
        
        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], s, 0, 0, [dzH], [], True)
        update_list_of_walls_and_volumes(Surfaces)
    #end of for loop

    #Connexion Droite
    Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], -xC1-longueur/2-R, h, 0, 0, 0, -1, math.pi/2, [dzC], [], True)
    update_list_of_walls_and_volumes(Surfaces)

    #Cylindre Vertical Droite
    Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[6], Surfaces[12], Surfaces[18]], 0, -h, 0, [dz], [], True)
    update_list_of_walls_and_volumes(Surfaces)

    print_GenMesh_log('fluid geometry done', sim_log)
    print_GenMesh_log('\n', sim_log)

    print_GenMesh_log('adding fluid physical surfaces ...', sim_log)
    gmsh.model.addPhysicalGroup(2, [face1, face2, face3, face4], name = "Inlet")
    gmsh.model.addPhysicalGroup(2, [Surfaces[0][1], Surfaces[6][1], Surfaces[12][1], Surfaces[18][1]], name = "Outlet")
    gmsh.model.addPhysicalGroup(2, list_of_walls, name = "fluid_to_solid")
    gmsh.model.addPhysicalGroup(3, list_of_volumes, name = "internalMesh")
    print_GenMesh_log('inlet, outlet, and walls patches created', sim_log)
    print_GenMesh_log('\n', sim_log)

    gmsh.model.geo.synchronize()
    print_GenMesh_log('generating fluid mesh ...', sim_log)
    gmsh.model.mesh.generate(3)
    print_GenMesh_log('fluid mesh generated', sim_log)
    print_GenMesh_log('\n', sim_log)
    print_GenMesh_log('writing fluid mesh file ...', sim_log)
    gmsh.write("fluid_mesh.msh2")
    print_GenMesh_log('fluid mesh file completed', sim_log)
    print_GenMesh_log('\n', sim_log)

    gmsh.finalize()
    ##############################################################

def generate_solid_mesh(plate_width, plate_height, plate_thickness, longueur, largeur, spacing, nombre_de_passes, perc_of_r, sim_log):
#plate_width, plate_height, longueur, largeur, spacing, nombre_de_passes, perc_of_r are used to draw the duct
#sim_log is used to continue filling the simulation's general log file
    #plate_width = plate width
    #plate_height : plate height
    #plate_thickness : plate thickness
    #longueur : duct width
    #largeur : duct thickness
    #spacing : distance between duct's borders and plate's borders on each side (left, right, and top)
    #nombre_de_passes : number of turns the duct makes inside the plate
    #perc_of_r : steepness of the connexion's curvature (the lower is perc_of_r, the stepper is the curvature)

    list_of_solid_walls.clear()

    space = plate_width - 2*spacing - 2*longueur
    height = plate_height - spacing - longueur

    R = perc_of_r * (longueur/2) /100

    xC1 = -space/2 #position de l'entrée
    h = height - longueur/2 - R #hauteur du canal d'entrée
    h_mid = h - (longueur/2+R) - longueur/2 -spacing #hauteur des canaux intermédiaires
    s = (space - 2*(longueur/2+R)*(2*nombre_de_passes-1)) / (2*nombre_de_passes-1) #espacements entre les canaux

    print_GenMesh_log('drawing solid geometry ...', sim_log)

    ##################Generation of solid Mesh####################
    gmsh.initialize()

    maillage = plate_thickness/15
    maillage_tube = maillage / 1 #finer mesh on tube-border

    #Bords de la Plaque
    point1 = gmsh.model.geo.add_point(-plate_width/2, 0, -plate_thickness/2, maillage)
    point2 = gmsh.model.geo.add_point(+plate_width/2, 0, -plate_thickness/2, maillage)
    point3 = gmsh.model.geo.add_point(+plate_width/2, 0, +plate_thickness/2, maillage)
    point4 = gmsh.model.geo.add_point(-plate_width/2, 0, +plate_thickness/2, maillage)

    point5 = gmsh.model.geo.add_point(-plate_width/2, plate_height, -plate_thickness/2, maillage)
    point6 = gmsh.model.geo.add_point(+plate_width/2, plate_height, -plate_thickness/2, maillage)
    point7 = gmsh.model.geo.add_point(+plate_width/2, plate_height, +plate_thickness/2, maillage)
    point8 = gmsh.model.geo.add_point(-plate_width/2, plate_height, +plate_thickness/2, maillage)

    #Bords du Tube Extérieur
    point10 = gmsh.model.geo.add_point(xC1-longueur/2, 0, -largeur/2, maillage_tube)
    point11 = gmsh.model.geo.add_point(xC1+longueur/2, 0, -largeur/2, maillage_tube)
    point12 = gmsh.model.geo.add_point(xC1+longueur/2, 0, +largeur/2, maillage_tube)
    point13 = gmsh.model.geo.add_point(xC1-longueur/2, 0, +largeur/2, maillage_tube)

    #Bords de la Plaque
    line1 = gmsh.model.geo.add_line(point1, point2)
    line2 = gmsh.model.geo.add_line(point2, point3)
    line3 = gmsh.model.geo.add_line(point3, point4)
    line4 = gmsh.model.geo.add_line(point4, point1)

    line5 = gmsh.model.geo.add_line(point5, point6)
    line6 = gmsh.model.geo.add_line(point6, point7)
    line7 = gmsh.model.geo.add_line(point7, point8)
    line8 = gmsh.model.geo.add_line(point8, point5)

    line9 = gmsh.model.geo.add_line(point1, point5)
    line10 = gmsh.model.geo.add_line(point2, point6)
    line11 = gmsh.model.geo.add_line(point3, point7)
    line12 = gmsh.model.geo.add_line(point4, point8)

    #Courbes du Tube Extérieur
    line13 = gmsh.model.geo.add_line(point10, point11)
    line14 = gmsh.model.geo.add_line(point11, point12)
    line15 = gmsh.model.geo.add_line(point12, point13)
    line16 = gmsh.model.geo.add_line(point13, point10)

    #Cylindre Vertical Gauche
    Surfaces = gmsh.model.geo.extrude([(1, line13), (1, line14), (1, line15), (1, line16)], 0, h, 0)
    update_solid_walls(Surfaces)

    #Connexion Gauche
    Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], xC1+longueur/2+R, h, 0, 0, 0, -1, math.pi/2)
    update_solid_walls(Surfaces)

    #Cylindre Horizontal
    Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], s, 0, 0)
    update_solid_walls(Surfaces)

    for i in range(nombre_de_passes - 1):
        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], xC1+(1+4*i)*(longueur/2+R)+(1+2*i)*s, h, 0, 0, 0, -1, math.pi/2)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], 0, -h_mid, 0)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], xC1+(3+4*i)*(longueur/2+R)+(1+2*i)*s, height-h_mid-(longueur/2+R), 0, 0, 0, -1, -math.pi/2)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], s, 0, 0)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], xC1+(3+4*i)*(longueur/2+R)+(2+2*i)*s, height-h_mid-(longueur/2+R), 0, 0, 0, -1, -math.pi/2)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], 0, h_mid, 0)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], xC1+(5+4*i)*(longueur/2+R)+(2+2*i)*s, h, 0, 0, 0, -1, math.pi/2)
        update_solid_walls(Surfaces)

        Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], s, 0, 0)
        update_solid_walls(Surfaces)
    #end of for loop

    #Connexion Droite
    Surfaces = gmsh.model.geo.revolve([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], -xC1-longueur/2-R, h, 0, 0, 0, -1, math.pi/2)

    update_solid_walls(Surfaces)

    #Cylindre Vertical Droite
    Surfaces = gmsh.model.geo.extrude([Surfaces[0], Surfaces[4], Surfaces[8], Surfaces[12]], 0, -h, 0)

    update_solid_walls(Surfaces)

    curve_loop1 = gmsh.model.geo.add_curve_loop([line1, line2, line3, line4])
    curve_loop11 = gmsh.model.geo.add_curve_loop([line13, line14, line15, line16])
    curve_loop12 = gmsh.model.geo.add_curve_loop([Surfaces[0][1], Surfaces[4][1], Surfaces[8][1], Surfaces[12][1]])
    curve_loop2 = gmsh.model.geo.add_curve_loop([line5, line6, line7, line8])
    curve_loop3 = gmsh.model.geo.add_curve_loop([line2, line11, -line10, -line6])
    curve_loop4 = gmsh.model.geo.add_curve_loop([line3, line12, -line11, -line7])
    curve_loop5 = gmsh.model.geo.add_curve_loop([line4, line9, -line8, -line12])
    curve_loop6 = gmsh.model.geo.add_curve_loop([line1, line10, -line5, -line9])

    face1 = gmsh.model.geo.add_plane_surface([curve_loop1, curve_loop11, curve_loop12])
    face2 = gmsh.model.geo.add_plane_surface([curve_loop2])
    face3 = gmsh.model.geo.add_plane_surface([curve_loop3])
    face4 = gmsh.model.geo.add_plane_surface([curve_loop4])
    face5 = gmsh.model.geo.add_plane_surface([curve_loop5])
    face6 = gmsh.model.geo.add_plane_surface([curve_loop6]) #heated face

    Surface_Loop_Faces = list(list_of_solid_walls)
    
    Surface_Loop_Faces.append(face1)
    Surface_Loop_Faces.append(face2)
    Surface_Loop_Faces.append(face3)
    Surface_Loop_Faces.append(face4)
    Surface_Loop_Faces.append(face5)
    Surface_Loop_Faces.append(face6)
    
    gmsh.model.geo.addSurfaceLoop(Surface_Loop_Faces, 1)
    gmsh.model.geo.addVolume([1], 1)

    print_GenMesh_log('solid geometry done', sim_log)
    print_GenMesh_log('\n', sim_log)

    print_GenMesh_log("adding solid physical surfaces ...", sim_log)
    gmsh.model.addPhysicalGroup(2, [face1, face2, face3, face4, face5], name = "Adiabatic_Walls")
    gmsh.model.addPhysicalGroup(2, [face6], name = "Heated_Face")
    gmsh.model.addPhysicalGroup(2, list_of_solid_walls, name = "solid_to_fluid")
    gmsh.model.addPhysicalGroup(3, [1], name = "internalMesh")
    print_GenMesh_log("adiabatic walls and heated face patches created", sim_log)
    print_GenMesh_log('\n', sim_log)

    gmsh.model.geo.synchronize()
    print_GenMesh_log('generating mesh', sim_log)
    gmsh.model.mesh.generate(3)
    print_GenMesh_log('mesh generated', sim_log)
    print_GenMesh_log('\n', sim_log)
    print_GenMesh_log('writing solid mesh file', sim_log)
    gmsh.write("solid_mesh.msh2")
    print_GenMesh_log('solid mesh file completed', sim_log)
    print_GenMesh_log('\n', sim_log)

    gmsh.finalize()
    ##############################################################


def Convert_Mesh_to_Foam_Format(simulation_name, sim_log):
    #Convert Meshes to Polymeshes
    print_GenMesh_log('converting fluid mesh from msh to foam format ...', sim_log)
    os.system('gmshToFoam fluid_mesh.msh2 -region fluid -case ' + simulation_name + ' > ' + simulation_name + '/fluid_mesh.log')
    print_GenMesh_log('convertion complete', sim_log)
    print_GenMesh_log('\n', sim_log)
    print_GenMesh_log('converting solid mesh from msh to foam format ...', sim_log)
    os.system('gmshToFoam solid_mesh.msh2 -region solid -case ' + simulation_name + ' > ' + simulation_name + '/solid_mesh.log')
    print_GenMesh_log('convertion complete', sim_log)
    os.system('rm -f fluid_mesh.msh2 solid_mesh.msh2')

    #adjusting fluid boundary file
    f = open(simulation_name + "/constant/fluid/polyMesh/boundary", "r")
    boundary = f.readlines()
    f.close()
    os.system('rm -f ' + simulation_name + '/constant/fluid/polyMesh/boundary')
    f = open(simulation_name + "/constant/fluid/polyMesh/boundary", "w")
    f.writelines(boundary[0:21])
    f.writelines(boundary[22:27])
    f.write('\t  type\t\t mappedWall;\n\t  inGroups\t\t 1(mappedPatch);\n')
    f.writelines(boundary[29:31])
    f.write('\t  sampleMode\t nearestPatchFace;\n\t  sampleRegion\t solid;\n\t  samplePatch\t solid_to_fluid;\n')
    f.writelines(boundary[31:35])
    f.writelines(boundary[36:43])
    f.close()

    #adjusting solid boundary file
    f = open(simulation_name + "/constant/solid/polyMesh/boundary", "r")
    boundary = f.readlines()
    f.close()
    os.system('rm -f ' + simulation_name + '/constant/solid/polyMesh/boundary')
    f = open(simulation_name + "/constant/solid/polyMesh/boundary", "w")
    f.writelines(boundary[0:20])
    f.write('\t  type\t\t mappedWall;\n\t  inGroups\t\t 1(mappedPatch);\n')
    f.writelines(boundary[22:24])
    f.write('\t  sampleMode\t nearestPatchFace;\n\t  sampleRegion\t fluid;\n\t  samplePatch\t fluid_to_solid;\n')
    f.writelines(boundary[24:28])
    f.writelines(boundary[29:35])
    f.writelines(boundary[36:43])
    f.close()

    print_GenMesh_log('\n', sim_log)
    
