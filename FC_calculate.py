#!/usr/bin/python

#USE: in the input file, set nstate higher by at least 2 than it would be necessary;
#calculate gradient, then modify mw_gradient_compact.txt to have only the
#states to be examined; then do a calculation with the original MINP file
#and write the results to "out.txt"; then call "cat <new_file>|python FC_calculate.py"

#LAST TWO STATES WON'T BE CALCULATED

import subprocess
import sys
import os

BOHR_TO_ANGSTROM = 0.529177249
STEP_SIZE = 0.6
NUMBER_OF_STEPS = 11


def call(command):
    #Calls an external command with shell=True.
    subprocess.call(command, shell=True)

def remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


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

def listifyGradient(string):
    #Makes a list containing the gradient out of a string
    list=string.split("====================\n")
    list.remove(list[0])
    for i in range(len(list)):
        list[i] = list[i].split("\n")
        #An empty string is added to list[i]; need to remove it:
        list[i] = list[i][0:-1]
    return(list)

def setMINP():
    #Creates new MINP file from the keywords and the modified geometry.
    with open("keywords","r") as keywords:
        keywords_contents=keywords.read()
    with open("geom","r") as geom:
        geom_contents=geom.read()
    remove("MINP")
    with open("MINP","a+") as MINP:
        MINP.write(keywords_contents)
        MINP.write(geom_contents)


def resetMINP():
    #Resets the MINP file.
    remove("MINP")
    os.rename("_MINP", "MINP")


def runMRCC(state,step):
    #Runs MRCC
    call("dmrcc > fc_out_"+str(state)+"_"+str(step)+".txt")
    output = open("fc_out_"+str(state)+"_"+str(step)+".txt","r")
    while True:
        line = output.readline()
        if not line:
            break
        if line.find("Error at the termination of mrcc.") != -1:
            raise Exception("""Something went wrong during calculation.
Please check fc_out_"""+str(state)+"""_"""+str(step)+""".txt and the input, and try again.""")
    output.close()



try:
    #save original MINP file in case something goes wrong
    call("cp MINP _MINP")

    #accepting the difference from the user (echo dr|python script.py)
    gradient_string = sys.stdin.read()
    gradient_list = listifyGradient(gradient_string)

    #getting the data
    input_data = get_input_data()
    calc = input_data["calc"]
    bohr_slash_unit = input_data["bohr_slash_unit"]
    nstate = input_data["nstate"]


    with open("original_geom","r") as file:
        original_string=file.read().rstrip("\n")
    geom_list=listifyCoordinates(original_string)

    #for every examined state
    for state in range(len(gradient_list)):
        #for every step
        for step in range(NUMBER_OF_STEPS):
            #reset geometry
            geom_list=listifyCoordinates(original_string)
            #for every atom:
            for i in range(len(geom_list)):
                #for every coordinate
                for j in range(1,4):
                    #add the coordinate and the corrseponding gradient * N_step * step_size
                    geom_list[i][j] = str(float(geom_list[i][j])-float(step)*STEP_SIZE*float(gradient_list[state][3*i+j-1])*bohr_slash_unit)
            new_string=delistifyCoordinates(geom_list)
            with open("geom","w") as file:
                file.write(new_string)
            setMINP()
            runMRCC(state,step)
            print("Calculation {0}, {1} done!".format(state,step))


    print("All done!")

finally:
    resetMINP()

    #deleting temporary files
    remove("original_geom")
    remove("geom")
    remove("keywords")
