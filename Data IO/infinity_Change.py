"""
The purpose of these commands was to fix an infinity symbol issue. 
I do not plan on making this a generic program, but something to reference to.
Created by Carlos A. Cartagena-Sanchez
"""
##################################################################################
""" You will have to adjust the max_Run_Number in line 29. 
    Also, you will have to change the filename variable on line 31-32."""
##################################################################################
import fileinput

def infinity_Change(filename, max_Value):

    with fileinput.FileInput(filename, inplace = True, backup = '.bak') as fi:
        for line in fi:
            print(line.replace(u"\u221E",str(max_Value)),end='')
            print(line.replace(u"\u221e",str(max_Value)),end='')
    
    return None

def Loop_Infinity(file_Path, pico_Number):
    ############################################################################################
    """ This particular script is ment to change the issue with the infinity symbol. 
        There is an issue and that is the max_Value. There are two types of data within
        the picoscope data files; magnetic and monitor values. Each will have a different
        max value. What one could do is write a script that separates the text file into
        two: magnetics and monitor values, then run two separate infinity scripts."""
    ############################################################################################
    max_Value = 1 #This is true for the magnetics, not the monitor values.
    max_Run_Number = 17 # for my current situation 04/23/2019
    for run in range(1,max_Run_Number + 1):
        filename = (file_Path + "04232019pico" + 
                    str(pico_Number) + "/20190423-0001 (" + str(run) + ").txt")
        infinity_Change(filename, max_Value)
    
    return None

file_Path = "04232019/"
##################################################################################
""" The line below can be looped for the four picoscopes."""
##################################################################################
for j in range(1,5):
    Loop_Infinity(file_Path, pico_Number = j)

