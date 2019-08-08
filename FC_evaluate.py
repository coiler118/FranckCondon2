#!/usr/bin/python

#USE: in the input file, set nstate higher by 2 than it would be necessary;
#calculate gradient, then modify mw_gradient_compact.txt to have only the
#states to be examined; then do a calculation with the original MINP file
#and write the results to "out.txt"; then call "cat <new_file>|python FC_calculate.py",
#then "cat <new_file>|python FC_evaluate.py"

#LAST TWO STATES WON'T BE CALCULATED

import subprocess
import sys
import os

BOHR_TO_ANGSTROM = 0.529177249
NUMBER_OF_STEPS = 11

#the name as which the total energy in the output can be found
NAME = {"CC2": "LR-CC2", "ADC(2)": "ADC(2)", "CCSD": "LR-CCSD", "SOS-ADC(2)": "ADC(2)", "SCS-ADC(2)": "ADC(2)"}

#DEFAULT SETTINGS
#maximum difference between two states to be recognized as each other
EQUAL_STATE_CRITERIUM = 4.0
#maximum energy difference between two MOs to be considered as close
MOS_CLOSE_CRITERIUM = 0.01
#maximum difference between the energy serial numbers of two states
#to be recognized as each other
MAX_ENERGY_DIFFERENCE = 3
#the weight of the difference between the transition overlap and 1.0
#in the difference between the state properties
TRANSITION_DIFFERENCE_FACTOR = 10.0
#maximum number of both natural and virtual MOs to be taken into account when
#examining the transitions
TRANSITION_MO_NUMBER = 20
#maximum number of dominant transitions to be printed
PRINTED_TRANSITIONS = 4

#FORCED_STATES format for example:
#"setting"
#"setting"
#1,3,5,6
#2,6,9,8
#...
#"1,3,5,6" means: while examining the energy of the molecule throughout
#the modification of its geometry according to the gradient of the 2nd (1+ always 1)
#examined state, force assign the state with the energy serial number 5 in step 3
#to the state with the energy serial number 6 in step 2 (3-1)
FORCED_STATES = []

#an option to modify the settings with a file
if os.path.exists("FCS"):
    FCS = open("FCS","r")
    charset="0123456789."
    while True:
        line=FCS.readline()
        if not line:
            break
        if line.find("EQUAL_STATE_CRITERIUM") != -1:
            num = []
            for char in line:
                if char in charset:
                    num.append(char)
            EQUAL_STATE_CRITERIUM = float("".join(num))
            set[0] = True
        elif line.find("MOS_CLOSE_CRITERIUM") != -1:
            num = []
            for char in line:
                if char in charset:
                    num.append(char)
            MOS_CLOSE_CRITERIUM = float("".join(num))
            set[1] = True
        elif line.find("MAX_ENERGY_DIFFERENCE") != -1:
            num = []
            for char in line:
                if char in charset:
                    num.append(char)
            MAX_ENERGY_DIFFERENCE = float("".join(num))

        elif line.find("TRANSITION_DIFFERENCE_FACTOR") != -1:
            num = []
            for char in line:
                if char in charset:
                    num.append(char)
            TRANSITION_DIFFERENCE_FACTOR = float("".join(num))

        elif line.find("TRANSITION_MO_NUMBER") != -1:
            num = []
            for char in line:
                if char in charset:
                    num.append(char)
            TRANSITION_MO_NUMBER = int("".join(num))

        elif line.find("PRINTED_TRANSITIONS") != -1:
            num = []
            for char in line:
                if char in charset:
                    num.append(char)
            PRINTED_TRANSITIONS = int("".join(num))


        elif line.find(",") != -1:
            line=line.strip()
            numbers_string = line.split(",")
            forced_state = tuple([int(i) for i in numbers_string])
            FORCED_STATES.append(forced_state)



class Transitions:
    #a class for the transitions
    def __init__(self, coeffs, Is, As, current_close_MOs):
        #coefficient of the transition
        self.coeffs = coeffs
        #the MO from which the transition happened
        self.Is = Is
        #the MO to which the transition happened
        self.As = As
        #the close MOs for the current examined state and step
        self.close_MOs = current_close_MOs

        #defining the length
        length_square = 0.0
        for coeff in self.coeffs:
            length_square += coeff*coeff
        self.length = length_square**0.5

        self.norm_coeffs = [coeff/self.length for coeff in self.coeffs]

        self.dominant_transition = []
        high_coeffs = self.coeffs
        for j in range(PRINTED_TRANSITIONS):
            highest_coeff = 0.0
            for i in range(len(high_coeffs)):
                if abs(high_coeffs[i]) > abs(highest_coeff):
                    highest_coeff = high_coeffs[i]
                    max_index = i
            if abs(highest_coeff) > 0.1:
                self.dominant_transition.append((self.Is[max_index], self.As[max_index]))
            high_coeffs[max_index] = 0.0

    def overlap(self,other):
        #calculate the overlap
        #if there are close MOs, calculate the possible products, and take the
        #largest one
        #doesn't take the signs into account
        overlap = 0.0
        #the number of examined virtual MOs
        virtual_MO_count = self.Is.count(self.Is[0])
        for i in range(len(self.coeffs)):
            possible_products = []
            possible_products.append(abs(self.norm_coeffs[i] * other.norm_coeffs[i]))

            for j in (-1,1):
                try:
                    if (self.As[i],self.As[i+j]) in self.close_MOs:
                        possible_products.append(abs(self.norm_coeffs[i]*other.norm_coeffs[i+j]))
                except IndexError:
                    pass

            for k in (-virtual_MO_count,virtual_MO_count):
                try:
                    if (self.Is[i], self.Is[i+k]) in self.close_MOs:
                        possible_products.append(abs(self.norm_coeffs[i] * other.norm_coeffs[i+k]))
                        for j in (-1,1):
                            try:
                                if (self.As[i],self.As[i+j]) in self.close_MOs:
                                    possible_products.append(abs(self.norm_coeffs[i]*other.norm_coeffs[i+j+k]))
                            except IndexError:
                                pass

                except IndexError:
                    pass

            overlap += max(possible_products)

        print(overlap)
        return overlap



def call(command):
    #Calls an external command with shell=True.
    subprocess.call(command, shell=True)

def remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


def listifyEnergies(energies):
    #Makes a list containing the energies out of a string.
    list=energies.split("\n")
    return list


def listifyGradient(string):
    #Makes a list containing the gradient out of a string
    list=string.split("====================\n")
    list.remove(list[0])
    for i in range(len(list)):
        list[i] = list[i].split("\n")
        #An empty string is added to list[i]; need to remove it:
        list[i] = list[i][0:-1]
    return(list)




def get_input_data():
    #gets calc, nstate, and unit from the MINP file
    with open("MINP","r") as file:
        file_list=file.readlines()
        bohr_slash_unit = BOHR_TO_ANGSTROM

        for count, line in enumerate(file_list):
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


    return {"calc": calc, "nstate": nstate, "bohr_slash_unit": bohr_slash_unit}


def getNumberOfOrbitals():
    #gets the number of the non-frozen natural orbitals
    #and the number of virtual orbitals
    with open("out.txt","r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            if line.find("Orbital energies [au]:") != -1:
                file.readline()
                natural = 0
                while True:
                    line = file.readline()
                    if line.find(".") == -1:
                        break
                    natural += 1
                virtual = 0
                while True:
                    line = file.readline()
                    if line.find(".") == -1:
                        break
                    virtual += 1

    with open("MINP","r") as file:
        number_of_H_atoms = 0
        while True:
            line = file.readline()
            if not line:
                break
            if line.find("geom") != -1:
                total_number_of_atoms = int(file.readline().strip())
                while True:
                    line = file.readline()
                    if not line:
                        break
                    if line.find("H") != -1:
                        number_of_H_atoms += 1
        number_of_heavy_atoms = total_number_of_atoms - number_of_H_atoms
    return [natural - number_of_heavy_atoms, virtual]


def getMOEnergyDifferences(file_handler):
    #Gets the differences between the energies of every two adjacent orbitals
    empty_row = 0
    while True:
        line = file_handler.readline()
        if empty_row == 2:
            break
        if line.find("Orbital energies [au]:") != -1:
            file_handler.readline()
            differences = []

            line = file_handler.readline()
            prev_MO_energy = []
            for i in range(16,27):
                prev_MO_energy.append(line[i])
            prev_MO_energy=float("".join(prev_MO_energy))

            while True:
                line = file_handler.readline()
                if line.find(".") == -1:
                    empty_row += 1
                    if empty_row == 2:
                        break
                else:
                    current_MO_energy = []
                    for i in range(16,27):
                        current_MO_energy.append(line[i])
                    current_MO_energy=float("".join(current_MO_energy))
                    difference = current_MO_energy - prev_MO_energy
                    differences.append(difference)
                    prev_MO_energy = current_MO_energy
    return differences

def getCloseMOs(every_MO_difference):
    #Determines which MOs are so close to each other in energy that their
    #energy order might have interchanged
    close_MOs = {}
    for state in range(len(every_MO_difference)):
        for step in range(NUMBER_OF_STEPS):
            close_MOs[(state,step)] = []
    #"state > step > MO difference" to "state > MO difference > step"
    ordered_list = []
    for state in range(len(every_MO_difference)):
        ordered_list.append([])
        for MO in range(len(every_MO_difference[0][0])):
            ordered_list[state].append([])
            for step in range(NUMBER_OF_STEPS):
                ordered_list[state][MO].append(every_MO_difference[state][step][MO])

    for state in range(len(ordered_list)):
        #examine the difference between two particular MOs for every step
        for MO in range(len(ordered_list[state])):
            if any([abs(difference) < MOS_CLOSE_CRITERIUM for difference in ordered_list[state][MO]]):
                close = True
            else:
                close = False

            if close:
                min_difference = 1.0
                min_index = -1
                for i in range(len(ordered_list[state][MO])):
                    if abs(ordered_list[state][MO][i]) < abs(min_difference):
                        min_difference = ordered_list[state][MO][i]
                        min_index = i
                #look for the step in which the two MOs are the closest, and
                #make it so that there is a possibility that the energy order
                #of the two MOs has interchanged in that and the two adjacent steps
                if min_index != 0 and min_index != NUMBER_OF_STEPS-1:
                    for d in range(-1,2):
                        close_MOs[(state,min_index+d)].append((MO+1,MO+2))
                        close_MOs[(state,min_index+d)].append((MO+2,MO+1))
                elif min_index == 0:
                    for d in range(0,2):
                        close_MOs[(state,min_index+d)].append((MO+1,MO+2))
                        close_MOs[(state,min_index+d)].append((MO+2,MO+1))
                elif min_index == NUMBER_OF_STEPS-1:
                    for d in range(-1,1):
                        close_MOs[(state,min_index+d)].append((MO+1,MO+2))
                        close_MOs[(state,min_index+d)].append((MO+2,MO+1))

    return close_MOs


def getStateProperties(file_handler,calc,energy_to_state,close_MOs,number_of_orbitals,relevant_state,step):
    #gets the energy serial number and the transition properties of a state
    properties = {}
    energy = 0
    while True:
        line = file_handler.readline()
        if not line:
            break
        if line.find("Total "+NAME[calc]+" energy") != -1:
            energy += 1
            state = energy_to_state[energy]
            properties[state] = [energy]
            for i in range(7):
                file_handler.readline()

            #getting the transition
            coefficients = []
            Is = []
            As = []
            count = -1
            while True:
                line = file_handler.readline()
                count += 1
                I_count = count//number_of_orbitals[1]  #number_of_orbitals[1] is the number of virtual orbitals
                A_count = count%number_of_orbitals[1]
                if line.find("============") != -1:
                    break
                #take the transition only if its MOs are close to the HOMO and the LUMO
                if I_count >= number_of_orbitals[0]-TRANSITION_MO_NUMBER and A_count < TRANSITION_MO_NUMBER:
                    #the coefficient of the transition
                    coefficient = []
                    for i in range(2,11):
                        coefficient.append(line[i])
                    coefficient = (float("".join(coefficient)))
                    coefficients.append(coefficient)

                    #the MO from which the transition happened
                    I = []
                    for i in range(14,17):
                        I.append(line[i])
                    I = int("".join(I))
                    Is.append(I)

                    #the MO to which the transition happened
                    A = []
                    for i in range(25,28):
                        A.append(line[i])
                    A = int("".join(A))
                    As.append(A)

            properties[state].append(Transitions(coefficients,Is,As,close_MOs[relevant_state,step]))


    return properties


def getStateEnergies(state,step,calc,prev_state_prop,nstate,close_MOs,number_of_orbitals):
    #Appends the total excited state energies to a file.
    #Determines which state in the previous step corresponds to
    #the energy serial number in the current step

    #get the energy serial number of the states and write them to a dictionary
    #(for the previous step)
    #energy serial number: state number
    energy_to_state_prev = {}
    for i in range(len(prev_state_prop)):
        energy_to_state_prev[prev_state_prop[i+1][0]] = i+1

    #get the current state properties
    #the state numbers will correspond to their energy serial numbers determined
    #by the previous step (for now)
    with open("fc_out_"+str(state)+"_"+str(step)+".txt","r") as output:
        current_state_prop = getStateProperties(output,calc,energy_to_state_prev,close_MOs,number_of_orbitals,state,step)

    #get the energy serial number of the states and write them to a dictionary
    #(for the current step)
    energy_to_state = {}
    for i in range(len(current_state_prop)):
        energy_to_state[current_state_prop[i+1][0]] = i+1

    new_state_prop = {}

    #get all excited states energies for the current step
    #add them to a list in the order of the energy, not the state
    energies_file = open("energies.txt","a+")
    output = open("fc_out_"+str(state)+"_"+str(step)+".txt","r")
    raw_energies = []
    while True:
        line = output.readline()
        if not line:
            break

        if line.find("Total "+NAME[calc]+" energy") != -1:
            line = line.strip()
            energy =[]
            colon_found = False
            for char in line:
                #append non-whitespace characters after the colon to energy_list
                if colon_found and char != " ":
                    energy.append(char)
                if char == ":":
                    colon_found = True
            energy = "".join(energy)
            raw_energies.append(energy)

    #print the state properties for the previous step and the current step
    print("\n")
    for i in range(1,nstate):
        print("{0}, {1}".format(prev_state_prop[i][0],prev_state_prop[i][1].dominant_transition))
        print("{0}, {1}".format(current_state_prop[i][0],current_state_prop[i][1].dominant_transition))
    print("\n")

    #examine the states in the order of their energies
    for current_energy in range(1,nstate):
        print("{0}, {1}".format(prev_state_prop[energy_to_state[current_energy]][0],prev_state_prop[energy_to_state[current_energy]][1].dominant_transition))
        print("{0}, {1}".format(current_state_prop[energy_to_state[current_energy]][0],current_state_prop[energy_to_state[current_energy]][1].dominant_transition))
        forced = False

        #forcing a state to be assigned to a state in the previous step
        for i in FORCED_STATES:
            if i[0] == state and i[1] == step and i[2] == current_energy:
                forced = True
                new_state_prop[energy_to_state[i[3]]] = current_state_prop[energy_to_state[current_energy]]
                #print the overlaps with the other states
                for other_energy in range(1,nstate):
                    overlap = current_state_prop[energy_to_state[current_energy]][1].overlap(prev_state_prop[energy_to_state[other_energy]][1])
                del overlap
                print("Forced {}".format(i[3]))

        if not forced:
            differences = []
            for other_energy in range(1,nstate):
                #defining the difference between the properties of two states
                difference = (1.0-current_state_prop[energy_to_state[current_energy]][1].overlap(prev_state_prop[energy_to_state[other_energy]][1]))*TRANSITION_DIFFERENCE_FACTOR+abs(current_energy-other_energy)
                differences.append(difference)
            #selecting the prev state with the lowest difference from the current one
            min_difference = 100.0
            for i in range(nstate-1):
                if differences[i] < min_difference and abs(i+1 - current_energy) <= MAX_ENERGY_DIFFERENCE:
                    min_difference = differences[i]
                    min_index = i
            if min_difference <= EQUAL_STATE_CRITERIUM:
                #assigning the state
                new_state_prop[energy_to_state[min_index+1]] =  current_state_prop[energy_to_state[current_energy]]
                print(min_index+1)
            #if the state has high energy, a new state may have appeared
            elif current_energy == nstate-1 or current_energy == nstate-2:
                new_state_prop[energy_to_state[current_energy]] = current_state_prop[energy_to_state[current_energy]]
                print("High energy state lost")
            #if no state has been found, stop the program
            else:
                raise Exception("Properties!")



    print("\n")
    for i in range(1,nstate):
        print("{0}, {1}".format(new_state_prop[i][0],new_state_prop[i][1].dominant_transition))
    print("\n")

    #creating the list of energies with the proper order (corresponding
    #to the state numbers)
    done_energies = []
    for i in range(nstate-3):
        done_energies.append(raw_energies[new_state_prop[i+1][0]-1])




    energies = "\n".join(done_energies)
    energies_file.write(energies+"\n")

    output.close()
    energies_file.close()

    #return the new state properties so it can be prev_state_prop for
    #the next step
    return new_state_prop

def orderEnergies(nstate,gradient_list):
    #Puts the energies in order: states to be examined > states > steps.


    remove("FC_energies.txt")

    #getting the energies as a string, and then as a list
    with open("energies.txt", "r") as energies_file:
        energies_string = energies_file.read()
    energies_list=listifyEnergies(energies_string)

    FC = open("FC_energies.txt", "a+")


    #for every examined state:
    for examined_state in range(len(gradient_list)):
        #Head for the examined state
        FC.write("\n====================\nExamined state ({}):\n".format(examined_state))

        #selecting the indices of energies belonging only to the examined state
        #len(energies_list) is the number of all calculated energies,
        #len(gradient_list) is the number of examined states
        #slicing the energies_list
        examined_state_energies = energies_list[len(energies_list)/len(gradient_list)*examined_state : len(energies_list)/len(gradient_list)*(examined_state+1)]
        for state in range(nstate-3):
            FC.write("\nState {}:\n".format(state+1))
            #indices belonging to the state
            state_energy_indices = [i for i in range(len(examined_state_energies)) if i%(nstate-3) == state]
            for step in range(NUMBER_OF_STEPS):
                #creating a one-element list containing the index of the energy belonging to
                #the state and the step
                step_index = [state_energy_indices[j] for j in range(len(state_energy_indices)) if j%NUMBER_OF_STEPS == step]

                energy = examined_state_energies[step_index[0]]
                FC.write(energy+"\n")
        FC.write("\n====================\n")

    FC.close()




try:
    print("""Settings used:
EQUAL_STATE_CRITERIUM = """+str(EQUAL_STATE_CRITERIUM)+"""
MOS_CLOSE_CRITERIUM = """+str(MOS_CLOSE_CRITERIUM)+"""
MAX_ENERGY_DIFFERENCE = """+str(MAX_ENERGY_DIFFERENCE)+"""
TRANSITION_DIFFERENCE_FACTOR = """+str(TRANSITION_DIFFERENCE_FACTOR)+"""
TRANSITION_MO_NUMBER = """+str(TRANSITION_MO_NUMBER)+"""
PRINTED_TRANSITIONS = """+str(PRINTED_TRANSITIONS)+"""
FORCED_STATES=""")
    for forced_state in FORCED_STATES:
        print(forced_state)
    #making a list out of the input gradient file
    gradient_string = sys.stdin.read()
    gradient_list = listifyGradient(gradient_string)

    input_data = get_input_data()
    calc = input_data["calc"]
    bohr_slash_unit = input_data["bohr_slash_unit"]
    nstate = input_data["nstate"]
    number_of_orbitals = getNumberOfOrbitals()

    #get the MO differences for every examined state and step
    every_MO_difference = []
    for state in range(len(gradient_list)):
        every_MO_difference.append([])
        for step in range(NUMBER_OF_STEPS):
            with open("fc_out_"+str(state)+"_"+str(step)+".txt","r") as output:
                differences = getMOEnergyDifferences(output)
                every_MO_difference[state].append(differences)

    #determine which MOs are close in which state and step
    close_MOs = getCloseMOs(every_MO_difference)

    #for every examined state
    for state in range(len(gradient_list)):

        #getting initial data
        initial_energies_to_state = {}
        for i in range(1,nstate):
            initial_energies_to_state[i] = i

        with open("out.txt", "r") as output:
            properties = getStateProperties(output,calc,initial_energies_to_state,close_MOs,number_of_orbitals,state,0)

        #for every step
        for step in range(NUMBER_OF_STEPS):
            properties = getStateEnergies(state,step,calc,properties,nstate,close_MOs,number_of_orbitals)
            print("Energies {0}, {1} done!".format(state,step))

    orderEnergies(nstate,gradient_list)

    print("All done!")

finally:

    #deleting temporary files
    remove("energies.txt")
