#!/usr/bin/python

#USE: call "echo <step size>|python GR_calculate.py"
import subprocess
import sys
import os

BOHR_TO_ANGSTROM = 0.529177249

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

def runMRCC(i,j,sign,calc):
    #Runs MRCC and appends the total excited state energies to a file.

    #the name as which the total energy in the output can be found
    name = {"CC2": "LR-CC2", "ADC(2)": "ADC(2)", "CCSD": "LR-CCSD"}

    #do the calculation
    call("dmrcc > gr_out_"+str(i)+"_"+str(j)+"_"+str(sign)+".txt")

    #raise error is something goes wrong
    with open("gr_out_"+str(i)+"_"+str(j)+"_"+str(sign)+".txt") as output:
        while True:
            line = output.readline()
            if not line:
                break
            if line.find("Error at the termination of mrcc.") != -1:
                raise Exception("""Something went wrong during calculation.
Please check gr_out_"""+str(i)+"""_"""+str(j)+"""_"""+str(sign)+""".txt and the input, and try again.""")




try:
    #save original MINP file in case something goes wrong
    call("cp MINP _MINP")

    #accepting the step size from the user (echo dr|script.py)
    dr = float(sys.stdin.read())

    #getting data from the input
    input_data = get_input_data()
    calc = input_data["calc"]
    bohr_slash_unit = input_data["bohr_slash_unit"]
    nstate = input_data["nstate"]


    with open("original_geom","r") as file:
        original_string=file.read().rstrip("\n")
    geom_list=listifyCoordinates(original_string)

    #for every atom:
    for i in range(len(geom_list)):
        #for every coordinate:
        for j in range(1,4):
            #modify geometry both negatively and positively:
            for sign in (-1,1):
                #reset geom_list
                geom_list=listifyCoordinates(original_string)
                #modify 1 coordinate, then make it a string
                geom_list[i][j]=str(float(geom_list[i][j])+(dr*float(sign)/2.0))
                new_string=delistifyCoordinates(geom_list)
                with open("geom","w") as file:
                    file.write(new_string)
                setMINP()
                runMRCC(i,j,sign,calc)
                print("Calculation {0}, {1}, {2} done!".format(i,j,sign))

    print("All done!")

finally:
    resetMINP()

    #deleting temporary files
    remove("original_geom")
    remove("geom")
    remove("keywords")
