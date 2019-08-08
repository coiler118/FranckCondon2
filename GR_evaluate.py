#!/usr/bin/python

#USE: call "echo <step size>|python GR_evaluate.py"
import subprocess
import sys
import os

BOHR_TO_ANGSTROM = 0.529177249
M_H = 1.00783
M_C = 12.00000
M_N = 14.00307
M_O = 15.99491
#the name as which the total energy in the output can be found
NAME = {"CC2": "LR-CC2", "ADC(2)": "ADC(2)", "CCSD": "LR-CCSD", "SOS-ADC(2)": "ADC(2)", "SCS-ADC(2)": "ADC(2)"}

def call(command):
    #Calls an external command with shell=True.
    subprocess.call(command, shell=True)

def remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


def listifyCoordinates(coordinates):
    #Makes a matrix containing the coordinates out of a string.
    list=coordinates.split("\n")
    for i in range(len(list)):
        list[i]=list[i].split()
    return list

def listifyEnergies(energies):
    #Makes a list containing the energies out of a string.
    list=energies.split("\n")
    return list

def delistifyCoordinates(list):
    #Makes a string out of the matrix containing the coordinates.
    newlist=[]
    #make a placeholder list
    for i in range(len(list)):
        newlist.append(i)
    #for every atom:
    for i in range(len(list)):
        newlist[i]="    ".join(list[i])
    coordinates="\n".join(newlist)
    return coordinates


def get_input_data():
    #gets calc, nstate, unit, the geometry and the keywords from the MINP file
    with open("MINP","r") as file:
        file_list=file.readlines()
        bohr_slash_unit = BOHR_TO_ANGSTROM
        remove("keywords")
        remove("original_geom")
        keywords = open("keywords","a+")
        geom = open("original_geom","a+")
        geom_reached = False
        #append the first 3 rows of the file to "keywords"
        for i in range(3):
            keywords.write(file_list[i])

        for count, line in enumerate(file_list):
            #if "geom" can be found in the line, set geom_reached to True and
            #get the next line
            if line.find("geom") != -1:
                geom_reached = True
                atoms_line = file_list[count+1]
            #append the 3rd following line to "keywords"
            #if "geom" has not been reached
            if not geom_reached:
                keywords.write(file_list[count+3])
            #append line to geom if it has coordinates
            elif line.strip() != "" and line.find("geom") == -1 and line != atoms_line:
                geom.write(line)

            #get calc (string following the equals sign)
            if line.find("calc") != -1:
                line = line.strip()
                calc_list = []
                equals_found = False
                for char in line:
                    if equals_found:
                        calc_list.append(char)
                    if char == "=":
                        equals_found = True
                calc = "".join(calc_list)

            #get nstate (number following the equals sign)
            elif line.find("nstate") != -1:
                line = line.strip()
                nstate_list = []
                equals_found = False
                for char in line:
                    if equals_found:
                        nstate_list.append(char)
                    if char == "=":
                        equals_found = True
                nstate = int("".join(nstate_list))

            #if the unit has been set to bohr, set bohr_slash_unit to 1.0
            elif line.find("bohr") != -1:
                bohr_slash_unit = 1.0

        keywords.close()
        geom.close()

    return {"calc": calc, "nstate": nstate, "bohr_slash_unit": bohr_slash_unit}



def getEnergies(i,j,sign,calc):
    #get the total energies and append them to a file
    with open("gr_out_"+str(i)+"_"+str(j)+"_"+str(sign)+".txt") as output:
        with open("energies.txt","a+") as energies:
            while True:
                line = output.readline()
                if not line:
                    break

                if line.find("Total "+NAME[calc]+" energy") != -1:
                    line = line.strip()
                    energy_list = []
                    colon_found = False
                    for char in line:
                        #append non-whitespace characters after the colon to energy_list
                        if colon_found and char != " ":
                            energy_list.append(char)
                        if char == ":":
                            colon_found = True
                    energy = "".join(energy_list)
                    energies.write(energy+"\n")

def calculateGradient(dr,bohr_slash_unit,nstate,geom_list):
    #Calculates the derivatives and puts them in order: states > atoms > coordinates

    #write the step size at the head of the file
    with open("gradient_elaborate.txt","w") as file:
        if bohr_slash_unit == 1.0:
            file.write("dr = "+str(dr)+" bohr\n")
        else:
            file.write("dr = "+str(dr)+" angstrom\n")

    with open("mw_gradient_elaborate.txt","w") as file:
        file.write("Mass-weighted\n")
        if bohr_slash_unit == 1.0:
            file.write("dr = "+str(dr)+" bohr\n")
        else:
            file.write("dr = "+str(dr)+" angstrom\n")

    remove("gradient_compact.txt")
    remove("mw_gradient_compact.txt")

    #getting the energies as a string, and then as a list
    with open("energies.txt", "r") as energies_file:
        energies_string = energies_file.read()
    energies_list=listifyEnergies(energies_string)

    elaborate = open("gradient_elaborate.txt", "a+")
    compact = open("gradient_compact.txt", "a+")
    #mass-weighted versions
    mw_elaborate = open("mw_gradient_elaborate.txt","a+")
    mw_compact = open("mw_gradient_compact.txt","a+")

    #for every excited state:
    for state in range(nstate-1):
        #Head for the state
        elaborate.write("\n======================================================================\nState {}:\n".format(state+1))
        mw_elaborate.write("\n======================================================================\nState {}:\n".format(state+1))
        compact.write("====================\n")
        mw_compact.write("====================\n")
        #selecting the indices of energies belonging only to the examined state
        state_energy_indices = [i for i in range(len(energies_list)) if i%(nstate-1)==state]
        #creating a list of the energy indices for both the positively and the negatively modified geometries
        positive_geom_indices = [state_energy_indices[j] for j in range(len(state_energy_indices)) if j%2 == 1]
        negative_geom_indices = [state_energy_indices[j] for j in range(len(state_energy_indices)) if j%2 == 0]
        #calculating the derivative
        for i in range(len(positive_geom_indices)):
            #write symbol of atom
            if i%3 == 0:
                symbol = geom_list[i/3][0]
                elaborate.write("\n"+symbol+"    ")
                mw_elaborate.write("\n"+symbol+"    ")

            derivative = ((float(energies_list[positive_geom_indices[i]]) - float(energies_list[negative_geom_indices[i]])) / dr * bohr_slash_unit)
            elaborate.write(str(derivative)+"    ")
            compact.write(str(derivative)+"\n")
            #atomic masses
            mass = {"H": M_H, "C": M_C, "N": M_N, "O": M_O}
            mw_elaborate.write(str(derivative*mass[symbol]**(-0.5))+"    ")
            mw_compact.write(str(derivative*mass[symbol]**(-0.5))+"\n")
        elaborate.write("\n======================================================================\n")
        mw_elaborate.write("\n======================================================================\n")

    elaborate.close()
    compact.close()
    mw_elaborate.close()
    mw_compact.close()





try:


    #accepting the step size from the user (echo dr|script.py)
    dr = float(sys.stdin.read())

    input_data = get_input_data()
    calc = input_data["calc"]
    bohr_slash_unit = input_data["bohr_slash_unit"]
    nstate = input_data["nstate"]


    with open("original_geom","r") as file:
        original_string=file.read().rstrip("\n")
    geom_list=listifyCoordinates(original_string)

    remove("energies.txt")
    #for every atom:
    for i in range(len(geom_list)):
        #for every coordinate
        for j in range(1,4):
            #modify geometry both negatively and positively
            for sign in (-1,1):
                getEnergies(i,j,sign,calc)
                print("Energies {0}, {1}, {2} done!".format(i,j,sign))


    calculateGradient(dr,bohr_slash_unit,nstate,geom_list)

    print("All done!")

finally:

    #deleting temporary files
    remove("original_geom")
    remove("geom")
    remove("keywords")
