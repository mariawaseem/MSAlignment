import sys

"""
CSE 190 
Homework 2: MS Alignment
Author: Maria Waseem
Start Date: 12/06/2023
"""


# Q1
# a
if sys.argv[1] == "q1a":
    aminoacidmasses =  {"A" : 71,
                        "R" : 156, 
                        "N" : 114, 
                        "D" : 115, 
                        "C" : 103, 
                        "E" : 129, 
                        "K" : 128, 
                        "G" : 57, 
                        "H" : 137, 
                        "L" : 113, 
                        "M" : 131, 
                        "F" : 147, 
                        "P" : 97, 
                        "S" : 87, 
                        "T" : 101, 
                        "W" : 186, 
                        "Y" : 163, 
                        "V" : 99 
                        }
    filename = sys.argv[2]
    sequence = sys.argv[3]
    modmass = int(sys.argv[4])
    modsequence = sequence[0:1] + "+" + str(modmass) + sequence[1:len(sequence)]
    file = open(filename, "r")
    spectrum = {1:0}
    newrow = []
    for line in file:
        newrow = line.split("\t")
        newrow[1] = int(newrow[1].strip())
        newrow[0] = int(newrow[0])
        spectrum[newrow[0]] = int(newrow[1])

    sequencemasses = []
    sequencemass = 0
    for i in range(len(sequence) + 1):
            sequencemass = 0
            for j in range(i):
                sequencemass += aminoacidmasses.get(sequence[j])
            sequencemasses.append(sequencemass)

    rows = len(spectrum) + 1
    cols = len(sequencemasses) + 1
    neginf = float('-inf')
    matrix = [[neginf] * cols for row in range(rows)]

    spectrumkeys = spectrum.keys()
    # p
    for col, mass in enumerate(sequencemasses, start = 1):
        matrix[0][col] = mass

    # s
    for row, key in enumerate(spectrumkeys, start = 1):
        matrix[row][0] = key


    modmatrix = [row[:] for row in matrix]

    # ----------------------------------------------------------------------------------------------------

    matrix[1][1] = 0
    modflag = False
    prevscore = neginf

    for i in range(2, rows):
        for j in range(2, cols):
            prevscores = []
            prevscores.append(0)
            for a in range (1, i):
                delta = (matrix[i][0] - matrix[i - a][0]) - (matrix[0][j] - matrix[0][j - 1])
                if delta > modmass:
                    break
                if (delta == 0):
                    if modflag:
                        if modmatrix[i - a][j - 1] != neginf:
                            prevscores.append(modmatrix[i - a][j - 1])
                    else:
                        if matrix[i - a][j - 1] != neginf:
                            prevscores.append(matrix[i - a][j - 1])
                elif ((delta == modmass) and not modflag):
                    if matrix[i - a][j - 1] != neginf:
                        modflag = True
                        prevscores.append(matrix[i - a][j - 1])

            if len(prevscores) == 0:
                continue

            prevscore = max(prevscores)
            if modflag:
                modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevscore
            else:
                matrix[i][j] = spectrum.get(matrix[i][0]) + prevscore

    last_column = [row[-1] for row in modmatrix[1:]]
    
    print(modsequence + " " + str(max(last_column)) + " 1")

# b
if sys.argv[1] == "q1b":
    aminoacidmasses =  {"A" : 71,
                        "R" : 156, 
                        "N" : 114, 
                        "D" : 115, 
                        "C" : 103, 
                        "E" : 129, 
                        "K" : 128, 
                        "G" : 57, 
                        "H" : 137, 
                        "L" : 113, 
                        "M" : 131, 
                        "F" : 147, 
                        "P" : 97, 
                        "S" : 87, 
                        "T" : 101, 
                        "W" : 186, 
                        "Y" : 163, 
                        "V" : 99 
                        }
    filename = sys.argv[2]
    sequence = sys.argv[3]
    modmass = int(sys.argv[4])
    file = open(filename, "r")
    spectrum = {1:0}
    newrow = []
    for line in file:
        newrow = line.split("\t")
        newrow[1] = int(newrow[1].strip())
        newrow[0] = int(newrow[0])
        spectrum[newrow[0]] = int(newrow[1])

    sequencemasses = []
    sequencemass = 0
    for i in range(len(sequence) + 1):
            sequencemass = 0
            for j in range(i):
                sequencemass += aminoacidmasses.get(sequence[j])
            sequencemasses.append(sequencemass)

    rows = len(spectrum) + 1
    cols = len(sequencemasses) + 1
    neginf = float('-inf')
    matrix = [[neginf] * cols for row in range(rows)]

    spectrumkeys = spectrum.keys()
    # p
    for col, mass in enumerate(sequencemasses, start = 1):
        matrix[0][col] = mass

    # s
    for row, key in enumerate(spectrumkeys, start = 1):
        matrix[row][0] = key


    modmatrix = [row[:] for row in matrix]

    # ----------------------------------------------------------------------------------------------------

    matrix[1][1] = 0
    prevbestscore = neginf
    maxmodpos = 0
    modflag = False

    for i in range(2, rows):
        for j in range(2, cols):
            prevbestscore = neginf
            for a in range(1, i):
                a = i - a
                delta = (matrix[i][0] - matrix[a][0]) - (matrix[0][j] - matrix[0][j - 1])
                if delta > modmass:
                    break
                if (delta == 0):
                    if modflag:
                        if modmatrix[a][j - 1] != neginf:
                            if (modmatrix[a][j - 1] > prevbestscore):
                                prevbestscore = modmatrix[a][j - 1]
                                modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevbestscore
                    else:
                        if matrix[a][j - 1] != neginf:
                            if (matrix[a][j - 1] > prevbestscore):
                                prevbestscore = matrix[a][j - 1]
                                matrix[i][j] = spectrum.get(matrix[i][0]) + prevbestscore
                elif (delta == modmass) and matrix[a][j - 1] > prevbestscore:
                    if matrix[a][j - 1] != neginf:
                        if (matrix[a][j - 1] > prevbestscore):
                            prevbestscore = matrix[a][j - 1]
                            modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevbestscore
                    maxmodpos = j
                    modflag = True

        if prevbestscore == neginf:
            continue

    last_column = [row[-1] for row in modmatrix[1:]]
    modsequence = sequence[0:maxmodpos - 1] + "+" + str(modmass) + sequence[maxmodpos - 1:len(sequence)]

    print(modsequence + " " + str(max(last_column)) + " " + str(maxmodpos - 1))

# Q2
# a
if sys.argv[1] == "q2a":
    aminoacidmasses =  {"A" : 71,
                        "R" : 156, 
                        "N" : 114, 
                        "D" : 115, 
                        "C" : 103, 
                        "E" : 129, 
                        "K" : 128, 
                        "G" : 57, 
                        "H" : 137, 
                        "L" : 113, 
                        "M" : 131, 
                        "F" : 147, 
                        "P" : 97, 
                        "S" : 87, 
                        "T" : 101, 
                        "W" : 186, 
                        "Y" : 163, 
                        "V" : 99 
                        }
    filename = sys.argv[2]
    sequence = sys.argv[3]
    file = open(filename, "r")
    spectrum = {1:0}
    newrow = []
    for line in file:
        newrow = line.split("\t")
        newrow[1] = int(newrow[1].strip())
        newrow[0] = int(newrow[0])
        spectrum[newrow[0]] = int(newrow[1])

    sequencemasses = []
    sequencemass = 0
    for i in range(len(sequence) + 1):
            sequencemass = 0
            for j in range(i):
                sequencemass += aminoacidmasses.get(sequence[j])
            sequencemasses.append(sequencemass)

    rows = len(spectrum) + 1
    cols = len(sequencemasses) + 1
    neginf = float('-inf')
    matrix = [[neginf] * cols for row in range(rows)]

    spectrumkeys = spectrum.keys()
    # p
    for col, mass in enumerate(sequencemasses, start = 1):
        matrix[0][col] = mass

    # s
    for row, key in enumerate(spectrumkeys, start = 1):
        matrix[row][0] = key


    modmatrix = [row[:] for row in matrix]
    matrix[1][1] = 0

    # ----------------------------------------------------------------------------------------------------

    highscore = 0
    highscoremod = 0

    for posmodmass in range(2, 51):
        rows = len(spectrum) + 1
        cols = len(sequencemasses) + 1
        neginf = float('-inf')
        prevbestscore = neginf
        maxmodpos = 0
        modflag = False
        
        for i in range(2, rows):
            for j in range(2, cols):
                prevbestscore = neginf
                for a in range(1, i):
                    a = i - a
                    delta = (matrix[i][0] - matrix[a][0]) - (matrix[0][j] - matrix[0][j - 1])
                    if delta > posmodmass:
                        break
                    if (delta == 0):
                        if modflag:
                            if modmatrix[a][j - 1] != neginf:
                                if (modmatrix[a][j - 1] > prevbestscore):
                                    prevbestscore = modmatrix[a][j - 1]
                                    modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevbestscore
                        else:
                            if matrix[a][j - 1] != neginf:
                                if (matrix[a][j - 1] > prevbestscore):
                                    prevbestscore = matrix[a][j - 1]
                                    matrix[i][j] = spectrum.get(matrix[i][0]) + prevbestscore
                    elif (delta == posmodmass) and matrix[a][j - 1] > prevbestscore:
                        if matrix[a][j - 1] != neginf:
                            if (matrix[a][j - 1] > prevbestscore):
                                prevbestscore = matrix[a][j - 1]
                                modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevbestscore
                        maxmodpos = j
                        modflag = True

            if prevbestscore == neginf:
                continue

        last_column = [row[-1] for row in modmatrix[1:]]
        if max(last_column) > highscore:
            highscore = max(last_column)
            highscoremod = posmodmass
            highscoremodmax = maxmodpos

    last_column = [row[-1] for row in modmatrix[1:]]
    modsequence = sequence[0:highscoremodmax - 1] + "+" + str(highscoremod) + sequence[highscoremodmax - 1:len(sequence)]

    print(modsequence + " " + str(highscore) + " " + str(highscoremodmax - 1))


# b
if sys.argv[1] == "q2b":
    aminoacidmasses =  {"A" : 71,
                        "R" : 156, 
                        "N" : 114, 
                        "D" : 115, 
                        "C" : 103, 
                        "E" : 129, 
                        "K" : 128, 
                        "G" : 57, 
                        "H" : 137, 
                        "L" : 113, 
                        "M" : 131, 
                        "F" : 147, 
                        "P" : 97, 
                        "S" : 87, 
                        "T" : 101, 
                        "W" : 186, 
                        "Y" : 163, 
                        "V" : 99 
                        }
    filename = sys.argv[2]
    sequence = sys.argv[3]
    file = open(filename, "r")
    spectrum = {1:0}
    newrow = []
    for line in file:
        newrow = line.split("\t")
        newrow[1] = int(newrow[1].strip())
        newrow[0] = int(newrow[0])
        spectrum[newrow[0]] = int(newrow[1])

    sequencemasses = []
    sequencemass = 0
    for i in range(len(sequence) + 1):
            sequencemass = 0
            for j in range(i):
                sequencemass += aminoacidmasses.get(sequence[j])
            sequencemasses.append(sequencemass)

    rows = len(spectrum) + 1
    cols = len(sequencemasses) + 1
    neginf = float('-inf')
    matrix = [[neginf] * cols for row in range(rows)]


    spectrumkeys = spectrum.keys()
    # p
    for col, mass in enumerate(sequencemasses, start = 1):
        matrix[0][col] = mass

    # s
    for row, key in enumerate(spectrumkeys, start = 1):
        matrix[row][0] = key


    modmatrix = [row[:] for row in matrix]
    matrix[1][1] = 0

    # ----------------------------------------------------------------------------------------------------

    highscore = 0
    highscoremod = 0

    for posmodmass in range(2, 51):
        rows = len(spectrum) + 1
        cols = len(sequencemasses) + 1
        neginf = float('-inf')

        prevbestscore = neginf
        maxmodpos = 0
        modflag = False
        
        for i in range(2, rows):
            for j in range(2, cols):
                prevbestscore = neginf
                for a in range(1, i):
                    a = i - a
                    delta = (matrix[i][0] - matrix[a][0]) - (matrix[0][j] - matrix[0][j - 1])
                    if delta > posmodmass:
                        break
                    if (delta == 0):
                        if modflag:
                            if modmatrix[a][j - 1] != neginf:
                                if (modmatrix[a][j - 1] > prevbestscore):
                                    prevbestscore = modmatrix[a][j - 1]
                                    modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevbestscore
                        else:
                            if matrix[a][j - 1] != neginf:
                                if (matrix[a][j - 1] > prevbestscore):
                                    prevbestscore = matrix[a][j - 1]
                                    matrix[i][j] = spectrum.get(matrix[i][0]) + prevbestscore
                    elif (delta== posmodmass) and matrix[a][j - 1] > prevbestscore:
                        if matrix[a][j - 1] != neginf:
                            if (matrix[a][j - 1] > prevbestscore):
                                prevbestscore = matrix[a][j - 1]
                                modmatrix[i][j] = spectrum.get(modmatrix[i][0]) + prevbestscore
                        maxmodpos = j
                        modflag = True

            if prevbestscore == neginf:
                continue

        last_column = [row[-1] for row in modmatrix[1:]]
        if max(last_column) > highscore:
            highscore = max(last_column)
            highscoremod = posmodmass
            highscoremodmax = maxmodpos

    last_column = [row[-1] for row in modmatrix[1:]]
    modsequence = sequence[0:highscoremodmax - 1] + "+" + str(highscoremod) + sequence[highscoremodmax - 1:len(sequence)]

    print(modsequence + " " + str(highscore) + " " + str(highscoremodmax - 1))

# Extra Credit
if sys.argv[1] == "ec":
    aminoacidmasses =  {"A" : 71,
                        "R" : 156, 
                        "N" : 114, 
                        "D" : 115, 
                        "C" : 103, 
                        "E" : 129, 
                        "K" : 128, 
                        "G" : 57, 
                        "H" : 137, 
                        "L" : 113, 
                        "M" : 131, 
                        "F" : 147, 
                        "P" : 97, 
                        "S" : 87, 
                        "T" : 101, 
                        "W" : 186, 
                        "Y" : 163, 
                        "V" : 99 
                        }
    
    filename = sys.argv[2]
    sequence = sys.argv[3]
    file = open(filename, "r")
    file.readline()
                
    spectrum = {1:0}
    newrow = []
    for line in file:
        newrow = line.split("\t")
        newrow[1] = int(newrow[1].strip())
        newrow[0] = int(newrow[0])
        spectrum[newrow[0]] = int(newrow[1])

    epicscore = 0
    epicscoremodpos = 0

    neginf = float('-inf')

    sequencemasses = []
    for i in range(len(sequence) + 1):
            sequencemass = 0
            for j in range(i):
                sequencemass += aminoacidmasses.get(sequence[j])
            sequencemasses.append(sequencemass)
    
    highscore = 0
    highscoremod = 0
    location = 0
    pred = []

    rows = len(spectrum) + 1
    cols = len(sequencemasses) + 1
    
    matrix = [[[0, [-1, -1, -1]] if (row == 1) else [neginf, [-1, -1, -1]] for col in range(cols)] for row in range(rows)]

    spectrumkeys = list(spectrum.keys())

    # p
    for col, mass in enumerate(sequencemasses, start = 1):
        matrix[0][col] = mass

    # s
    for row, key in enumerate(spectrumkeys, start = 1):
        matrix[row][0] = key

    modmatrix = [row[:] for row in matrix]
    
    for posmodmass in range(2, 51):
        maxmodpos = 0
        modflag = False
        
        for i in range(2, rows):
            for j in range(2, cols):
                prevbestscore = neginf
                for a in range(1, i):
                    a = i - a
                    delta = (matrix[i][0] - matrix[a][0]) - (matrix[0][j] - matrix[0][j - 1])
                    if delta > posmodmass:
                        break
                    if (delta == posmodmass) and matrix[a][j - 1][0] > prevbestscore:
                        if matrix[a][j - 1] != neginf:
                            if (matrix[a][j - 1][0] > prevbestscore):
                                prevbestscore = matrix[a][j - 1][0]
                                modmatrix[i][j] = [spectrum.get(modmatrix[i][0]) + prevbestscore, [a, j-1, 0]]
                                if modmatrix[i][j][0] > highscore:
                                    highscore = modmatrix[i][j][0]
                                    highscoremod = posmodmass
                                    highscoremodmax = j
                                    location = j
                                    pred = modmatrix[i][j][1]
                        modflag = True
                    elif (delta == 0):
                        if modflag:
                            if modmatrix[a][j - 1][0] != neginf:
                                if (modmatrix[a][j - 1][0] > prevbestscore):
                                    prevbestscore = modmatrix[a][j - 1][0]
                                    modmatrix[i][j] = [spectrum.get(modmatrix[i][0]) + prevbestscore, [a, j-1, 1]]
                                    if modmatrix[i][j][0] > highscore:
                                        highscore = modmatrix[i][j][0]
                                        highscoremod = posmodmass
                                        location = j
                                        pred = modmatrix[i][j][1]
                        else:
                            if matrix[a][j - 1][0] != neginf:
                                if (matrix[a][j - 1][0] > prevbestscore):
                                    prevbestscore = matrix[a][j - 1][0]
                                    matrix[i][j] = [spectrum.get(matrix[i][0]) + prevbestscore, [a, j-1, 0]]
                                    if matrix[i][j][0] > highscore:
                                        highscore = matrix[i][j][0]
                                        highscoremod = posmodmass
                                        location = j
                                        pred = matrix[i][j][1]

            if prevbestscore == neginf:
                continue

    cell = modmatrix[pred[0]][pred[1]]
    positionofmod = 0
    count = 0
    while cell[1][2] == 1:
        cell = modmatrix[cell[1][0]][cell[1][1]]
        positionofmod += 1
        count += 1
    while cell[1][2] == 0:
        cell = matrix[cell[1][0]][cell[1][1]]
        count += 1

    startaa = location - count - 2
    epicscoremodpos = location - positionofmod - 1
    
    modsequence = sequence[startaa:epicscoremodpos - 1] + "+" + str(highscoremod) + sequence[epicscoremodpos - 1: location - 1]

    print(modsequence + " " + str(highscore) + " " + str(count - positionofmod))
    