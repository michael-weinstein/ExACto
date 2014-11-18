#!/usr/bin/env python3
'''
Exome aggregation consortium annotation tool, Version 0.3
This program is designed to use the ExAC VCF data from the broad institute (downloadable from the link below).
ftp://ftp.broadinstitute.org/pub/ExAC_release
You will first need to run the library generation function using the commandline argument -s filename.vcf
After that, you will be able to annotate your variants using the commandline argument -f filename.txt
This version has not yet been tested against an ANNOVAR output, so it will probably encounter errors there.
Copyright 2014, Michael Weinstein, Cohn Lab, UCLA.  This is free to use for noncommercial entities, but please
let me know that you are using if I don't already know.  mweinste@ucla.edu
To do:  Check items going in the output string to see if they need quotes and quote them.
'''

def checkargs():  #subroutine for validating commandline arguments
    import argparse #loads the required library for reading the commandline
    import os  #imports the library we will need to check if the file exists
    file = False
    job = False
    parser = argparse.ArgumentParser()
    parser.add_argument ("-f", "--file", help = "Specify the desired file to annotate for submission.")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("-s", "--split", help = "Specify a VCF to split into a subVCF library for use as a reference.")
    args = parser.parse_args()  #puts the arguments into the args object
    if not (args.file or args.split):  #if the args.file value is null, give an error message and quit the program
        usage("No file specified for annotation or splitting.") 
        quit()
    elif (args.file and args.split):
        quit('You must specify either a file to split or a file to annotate.  You cannot specify both on a single run.  Please try again with one at a time.\nGoodbye.')
    elif args.file:
        file = args.file
        job = 'annotate'
        if not os.path.isfile(file):
            quit('Specified file for annotation does not exist')
        if os.path.isfile(file + 'ExACtoed.txt'):
            quit('Specified file for annotation already has an existing annotation from this program.  Please move or delete it, then run again.')
    elif args.split:
        file = args.split
        job = 'librarysplit'
        if not os.path.isfile(file):
            quit('Specified VCF for the library split does not exist.')
        if os.path.exists('subvcfs'):
            quit('Library directory already exists.  Please remove the old one before creating a new one')
    return (file, job)

def usage(sin):  #This subroutine prints directions
    print ('Error: ' + sin)
    print ('This program will annotate a SeattleSeq or ANNOVAR (not yet tested) output and annotate the variants based upon the Exome Aggregation Consortium\'s VCF data.  Before your first annotation run, you will need to download the current ExAC VCF and run the library generation on it using the commandline argument -s ExACVCF.vcf.  Once that is complete, you can annotate your outputs using the commandline argument -f yourfile.txt.')

def directoryhash(directory):
    import os  #the library needed to look at files and directories
    import hashlib #the library needed to create an MD5 hash of the directory
    try:
        allfiles = os.listdir(directory)  #gets a list of all the files in the directory
        filestring = ''.join(allfiles)  #turns the list of files into a single string
        filestringhash = hashlib.md5(rawbytes(filestring)) #turns the string of filenames into an MD5 hash
    except:  #if something gives an error here
        return False #return a value of false instead of returning the hash string to let the program know it did not work
    return filestringhash

def rawbytes(stringin):  #subroutine to take a string and make it into UTF-8 bytes, often for sending via socket connection
    return bytes(stringin, 'utf-8')

def checkintegrity(directory):
    try:
        hashfile = open(directory + '/hashsum', 'r')
        recordedhash = (hashfile.readline()).strip('\n')
        hashfile.close()
        if recordedhash and recordedhash == (directoryhash(directory)).hexdigest():
            return True
    except:
        return False
    else:
        return False
        
def createintegrity(directory):
    try:
        hashfile = open(directory + '/hashsum', 'w')  #opens a file to save the recorded hash at time of creation
        checksum = directoryhash(directory).hexdigest() #creates a hash of the directory filenames and saves it as checksum
        if not checksum:  #if we did not get back a checksum value (because the subroutine failed)
            return False #return false to let the program know it did not work
        hashfile.write(checksum)  #otherwise, we write the checksum to a file
        hashfile.close() #and then close it
    except:  #if the process fails in some unexpected way
        return False  #return false to let the rest of the program know that
    return True  #otherwise return true to indicate success

def librarysplit(filename):  #subroutine for generating the subvcf library to speed searches
    import re  #the library we need to use regular expressions
    import os  #the library we need to look for files and directories
    if os.path.exists('subvcfs'):  #check if there is already a subvcfs directory, if so, we do not want to overwrite it
        usage('Library directory already exists.  Please remove the old one before creating a new one')
        quit()
    try:
        os.makedirs('subvcfs')  #actually creates the new directory
    except PermissionError: #if it is trying to create a directory and lacks the permission, we give an error message and quit
        print ('\nPermission to create the library directory was denied.  Please run this program with appropriate permissions or create it using your file manager.')
        quit()
    except:  #if we get some other error, give a message and quit
        print ('\nUnable to make the library directory.')
        quit()    
    lastchromosome = 0  #initializing a bunch of values here
    linecount = 0
    libraryfile = False
    lastpositionblock = 0
    try:  #see if we can open the exac VCF and read a line from it
        exac = open(filename, 'r')
        line = exac.readline()
    except: #otherwise, if we get an error, we display an error message and quit
        quit('Error reading from ExAC VCF')
    while line:  #while line is not empty (indiciating we have hit the end of the file)
        linecount += 1  #increment the linecount variable
        if re.match('^\#', line):  #check if the line starts with a hashtag, indicating that it is not a data line
            try:  #try to read the next line and then restart the loop
                line = exac.readline()
                continue
            except: #if the reading fails, display an error message
                quit('Error reading from ExAC VCF')
        linearray = line.split('\t')  #if the line did not start with a hashtag (indicating it is a data line), we will split it at every tab
        chromosome = linearray[0] #finds the chromosome
        position = linearray[1] #finds the position
        positionblock = str(int(position) - (int(position) % 100000)) #rounds down the position value to the nearest 100,000 bases to know which library file to use
        if positionblock != lastpositionblock or chromosome != lastchromosome:  #if either the chromosome or position block for this variant is different from the last
            if libraryfile:  #if we already have an open subvcf that we were writing to, we close it
                libraryfile.close()
            libraryfile = createlibraryfile(chromosome, positionblock) #then we create a new library file for the current chromosome and position block combination
        try:
            libraryfile.write(line)  #writes the data line to the subvcf.  By this point, we either know that we are writing to the same file as the last iteration, or have already opened the new one for writing
        except:
            quit('Error writing new library file.')  #if something goes wrong with writing, we show an error message
        try:
            line = exac.readline() #then we try to read the next line from the file
        except:
            quit('Error reading from ExAC VCF') #and display an error message if something goes wrong
        lastpositionblock = positionblock  #remember which chromosome and position block we were on during this iteration so we can compare with the next
        lastchromosome = chromosome
        print('Processed ' + str(linecount) + ' lines.', end = '\r')  #and tell the user what the progress is
    print('\n')
    if not createintegrity('subvcfs'):  #if we were not able to create the hashsum file to monitor the integrity of the library directory
        return False #returns a value of false to let the rest of the program know
    return True #otherwise returns a value of true, indicating success

def createlibraryfile(chromosome, positionblock): #subroutine for generating subvcf files for the library
    libraryfilename = str('subvcfs/' + chromosome + 'c' + positionblock + '.subvcf')  #uses the chromosome and position block to generate the actual filename
    try:  #try to create a new file using that filename
        libraryfile = open(libraryfilename, 'w') #needs a try/except block
    except: #if that fails, display an error message and quit
        quit('Error creating new library file.')
    return libraryfile #returns the newly created file handle to the part of the program that called it

def yesanswer(question):  #asks the question passed in and returns True if the answer is yes, False if the answer is no, and keeps the user in a loop until one of those is given.  Also useful for walking students through basic logical python functions
    answer = False  #initializes the answer variable to false.  Not absolutely necessary, since it should be undefined at this point and test to false, but explicit is always better than implicit
    while not answer:  #enters the loop and stays in it until answer is equal to True
        print (question + ' (Y/N)')  #Asks the question contained in the argument passed into this subroutine
        answer = input('>>') #sets answer equal to some value input by the user
        if str(answer) == 'y' or str(answer) == 'Y':  #checks if the answer is a valid yes answer
            return True  #sends back a value of True because of the yes answer
        elif str(answer) == 'n' or str(answer) == 'N': #checks to see if the answer is a valid form of no
            return False  #sends back a value of False because it was not a yes answer
        else: #if the answer is not a value indicating a yes or no
            print ('Invalid response.')
            answer = False #set ansewr to false so the loop will continue until a satisfactory answer is given

def quotedsplit(input, delimiter = ','):  #subroutine for splitting CSV lines (that will often have quotes around them and be separated by commas)  This line can take an optional argument of something other than a comma as a delimiter
    import re
    brokenlist = input.split(delimiter)  #this will split the input at every delimiter, quoted or not
    outputlist = [] #initialize a blank list for the output
    outputindex = -1  #initializing the index counter for the output list to -1 (because lists start at 0 for the first element)
    for item in brokenlist:
        if re.match('^".*[^"]$', outputlist[outputindex]):  #if the item in the current cell in the output list starts with quotes, but does not end with them
            outputlist[outputindex] = outputlist[outputindex] + delimiter + item  #append the next item from the broken, separated by a delimiter (fixes the broken entries)
        else:
            outputlist[outputindex] = outputlist[outputindex].strip('"') #removes the quotes from the entry in the list (they only get in the way at this point)  we will have to return them before writing them to a file
            outputindex += 1  #if the item in the output list looks complete
            outputlist.append(item) #create a new cell in the list containing the next item from the broken list (we will check to see if it's broken on the next iteration)
    return outputlist  #returns the list to the program for use
    
def getkeycolumns(headers):
    import re
    chromosomecolumn = False
    positioncolumn = False
    for i in range (0, len(headers)):
        if re.search('^chr$', headers[i], re.IGNORECASE) or re.search('^chromosome$', headers[i], re.IGNORECASE):
            chromosomecolumn = i
        elif re.search('^start$', headers[i], re.IGNORECASE) or re.search('^position$', headers[i], re.IGNORECASE):
            positioncolumn = i
        elif re.search('^ref$', headers[i], re.IGNORECASE) or re.search('^referencebase$', headers[i], re.IGNORECASE):
            reference = i
        elif re.search('^obs$', headers[i], re.IGNORECASE) or re.search('^samplealleles$', headers[i], re.IGNORECASE):
            observed = i
    if not chromosomecolumn and positioncolumn and reference and observed:
        return False
    return (chromosomecolumn, positioncolumn, reference, observed)

def createzerohash(populations):
    multiallele = False
    datapoints = ['AC','AN','Hom']
    fillins = {'AC':0,'AN':1,'Hom':0}
    frequencyhash = {}
    for point in datapoints:    
        for population in populations:
            try:
                frequencyhash[point][population] = fillins[point]
            except KeyError:
                frequencyhash[point]={}
                frequencyhash[point][population] = fillins[point]
    return frequencyhash

def getfrequencyhash(line, populations, position):
    import re
    multiallele = False
    populationset = ''
    frequencyhash = {}
    for population in populations:
        if populationset:
            populationset += '|'
        populationset += population
    regex = re.compile('(AC|AN|Hom)_(' + populationset + ')=(.+?);')
    dataset = re.findall(regex, line)
    dataset = [list(entry) for entry in dataset]  #turns a tupple into a list so that it can be modified properly later (this is not an issue unless we have a multi-allele variant)
    for datum in dataset:
        if ',' in datum[2]:
            multiallele = True
            break
    if multiallele:
        for i in range (0,len(dataset)):
            if ',' in dataset[i][2]:
                dataset[i][2] = dataset[i][2].split(',')
    for datum in dataset:
        if multiallele and datum[0] != 'AN':
            try:
                frequencyhash[datum[0]][datum[1]] = datum[2][position]
            except KeyError:
                frequencyhash[datum[0]]={}
                frequencyhash[datum[0]][datum[1]] = datum[2][position]
        else:
            try:
                frequencyhash[datum[0]][datum[1]] = datum[2]
            except KeyError:
                frequencyhash[datum[0]]={}
                frequencyhash[datum[0]][datum[1]] = datum[2]
    return frequencyhash

def annotate(file):
    import os
    import re
    if not os.path.isdir('subvcfs'):
        usage('No subvcf library detected in this directory.  Please be sure the script is running from the same directory that contains the library.')
        quit()
    fileformat = 'tdt'
    delimiter = '\t'
    lastlibraryfile = False
    populations = ['AFR','AMR','EAS','FIN','NFE','SAS']
    datapoints = ['AF','HOMOF']
    if not checkintegrity('subvcfs'):
        quit('Check of subvcf library directory not passed.  Please be sure no files have been added, removed, or renamed in the library director')
    inputfile = open(file, 'r')
    try:
        header = inputfile.readline().strip('\r\n\t')
    except:
        quit('Error reading header from input file.')
    if not '\t' in header:
        if not yesanswer('This does not appear to be a tab-delimited file.  Is it comma-delimited?'):
            quit('Please talk to Mike about getting this script to take your file format.')
        else:
            fileformat = 'csv'
    if fileformat == 'tdt':
        headerlist = header.split('\t')
        delimiter = '\t'
    if fileformat == 'csv':
        headerlist = quotedsplit(header)
        delimiter = ','
    headercolumns = len(headerlist)
    keycolumns = getkeycolumns(headerlist)
    if not keycolumns:
        quit('Unable to find columns containing chromosome and/or position data.')
    chromosomecolumn = keycolumns[0]
    positioncolumn = keycolumns[1]
    referencecolumn = keycolumns[2]
    observedcolumn = keycolumns[3]
    outputfilename = file + '.ExACtoed.txt'
    outputfile = open(outputfilename, 'w', 1)
    try:
        line = inputfile.readline().strip('\r\n\t')
    except:
        quit('Error reading first data line from input file.')
    linenumber = 0
    outputgroups = []
    for population in populations:
        outputgroups.append(population)
    outputgroups.append('max')
    outputheaderlist = []
    for header in headerlist:
        outputheaderlist.append(header)
    for allele in range(0,2):
        for point in datapoints:
            for group in outputgroups:
                outputheaderlist.append('A' + str(int(allele) + 1) + '_' + point + '_' + group)
    outputheaderlist.append('F_rarest_allele')
    outputheaderlist.append('Combo_max')
    outputheaderstring = ''
    for item in outputheaderlist:
        if outputheaderstring:
            outputheaderstring += delimiter
        outputheaderstring += item
    outputheaderstring += '\n'
    try:
        outputfile.write(outputheaderstring)
    except:
        quit('Error writing header to output file.')
    lastchromosome = False
    lastposition = False
    lastreference = False
    lastobserved = False
    newlines = 0
    repeatlines = 0
    while line:
        line = line.strip('\r\n\t')
        linenumber += 1
        if line[0] == '#':
            print('Skipped line ' + str(linenumber) + ' as it appears to be a comment.')
            try:
                outputfile.write(line)
            except:
                quit('Error writing to output file.')
            line = inputfile.readline()
            continue
        if linenumber == -1:  #this if/then statement here entirely for debugging purposes and nothing else
            pass
        print ('Processing line ' + str(linenumber), end = '\r')
        if fileformat == 'tdt':
            linearray = line.split('\t')
        elif fileformat == 'csv':
            linearray = quotedsplit(line)
        chromosome = linearray[chromosomecolumn]
        position = linearray[positioncolumn]
        reference = linearray[referencecolumn]
        observed = linearray[observedcolumn]
        if not (chromosome == lastchromosome and position == lastposition and reference == lastreference and observed == lastobserved): #In an output listing several possible transcripts for each variant, this recycling of already-calculated values cuts more than 50% of the time required.
            newlines += 1
            frequencyhash = {}
            refmismatch = False
            extraalleles = False
            observedarray = observed.split('/')
            currentlibraryfile = str(chromosome) + 'c' + str(int(position) - (int(position) % 100000)) + '.subvcf'
            if os.path.isfile('subvcfs/' + currentlibraryfile):
                founddata = True
            else:
                founddata = False
            if founddata:
                if currentlibraryfile != lastlibraryfile:
                    try:
                        libraryfile = open('subvcfs/' + currentlibraryfile, 'r')
                        library = libraryfile.read() #slurps the whole file into the library string
                        libraryfile.close()
                    except:
                        quit('Error opening subvcf library file.')
                refline = re.search('^(' + str(chromosome) + '\t' + str(position) + '\t.*?)$', library, re.MULTILINE)
                if not refline:
                    founddata = False
                else:
                    refline = refline.group(0)
                    reflinearray = refline.split('\t')
                if founddata:
                    exacalts = reflinearray[4].split(',')
                    exacalthash = {}
                    for i in range (0, len(exacalts)):
                        exacalthash[exacalts[i]]= i
                    if refline and (reflinearray[3] != reference):  #if this doesn't analyze lazy, we need two if statements
                        print('Reference base mismatch starting on line ' + str(linenumber) + '.')
                        refmismatch = True
                    else:
                        try:
                            observedarray.remove(reference)
                            singlenonref = True
                        except:
                            singlenonref = False
                            if len(observedarray) == 1:
                                homozygousrare = True
                            elif len(observedarray) > 1 and observedarray[0] == observedarray[1]: #if this doesn't lazy evaluate, it needs to be nested
                                homozygousrare = True
                                observedarray.remove(observedarray[0])
                            else:
                                homozygousrare = False
                        if len(observedarray) in range(1,3):
                            extraalleles = False
                            for allele in (observedarray):
                                try:
                                    frequencyhash[allele] = getfrequencyhash(reflinearray[7], populations, exacalthash[allele])
                                except KeyError:
                                    frequencyhash[allele] = createzerohash(populations)
                        else:
                            print('Possible error on line ' + str(linenumber) + ', subject appears to have more than 2 alleles.')
                            extraalleles = True
                        if not extraalleles:
                            i = 0
                            namelist = []
                            namehash = {}
                            valuearray = []
                            for allele in observedarray:
                                for point in datapoints:
                                    for population in populations:
                                        try:
                                            namehash[allele][point][population] = i
                                        except KeyError:
                                            try:
                                                namehash[allele][point] = {}
                                                namehash[allele][point][population] = i
                                            except KeyError:
                                                try:
                                                    namehash[allele] = {}
                                                    namehash[allele][point] = {}
                                                    namehash[allele][point][population] = i
                                                except:
                                                    print('This is what you get for not making the data structure simpler')
                                        namelist.append([allele, point, population])
                                        if point == 'AF':
                                            try:
                                                entry = int(frequencyhash[allele]['AC'][population])/int(frequencyhash[allele]['AN'][population])
                                            except KeyError:
                                                entry = 0
                                            except ZeroDivisionError:
                                                entry = 'NA'
                                        elif point == 'HOMOF':
                                            if chromosome not in ('X','x','Y','y','MT','Mt','mt'):
                                                try:
                                                    entry = int(frequencyhash[allele]['Hom'][population])/(int(frequencyhash[allele]['AN'][population])/2)
                                                except KeyError:
                                                    entry = 0
                                                except ZeroDivisionError:
                                                    entry = 'NA'
                                            else:
                                                entry = 'NA'
                                        valuearray.append(entry)
                                        i += 1
                                    maximum = 0
                                    namehash[allele][point]['max'] = i
                                    namelist.append([allele, point, 'max'])
                                    for population in populations:
                                        if valuearray[namehash[allele][point][population]] != 'NA':
                                            if maximum < valuearray[namehash[allele][point][population]]:
                                                maximum = valuearray[namehash[allele][point][population]]
                                    valuearray.append(maximum)
                                    i += 1
                            if not singlenonref and not homozygousrare:
                                allele1max = valuearray[namehash[observedarray[0]][point]['max']]
                                allele2max = valuearray[namehash[observedarray[1]][point]['max']]
                                rarestallele = min(int(allele1max),int(allele2max))
                                if allele1max == 0:
                                    allele1max = 1/65000
                                if allele2max == 0:
                                    allele2max = 1/65000
                                combomax = allele1max * allele2max
                            elif singlenonref:
                                rarestallele = valuearray[namehash[observedarray[0]][point]['max']]
                                combomax = 'NA'
                            elif homozygousrare:
                                rarestallele = valuearray[namehash[observedarray[0]][point]['max']]
                                combomax = (1/65000)**2
                            datastring = ''
                            for i in range(0,headercolumns):
                                if datastring:
                                    datastring += delimiter
                                if delimiter in linearray[i]:
                                    linearray[i] = '"' + linearray[i] + '"'
                                datastring += linearray[i]
                            for value in valuearray:
                                datastring += delimiter + str(value)
                            if len(observedarray) == 1:
                                for value in valuearray:
                                    datastring += delimiter + 'NA'
                                datastring += delimiter + str(valuearray[namehash[observedarray[0]][point]['max']]) + delimiter + str(valuearray[namehash[observedarray[0]][point]['max']]**2)
                            for i in range (headercolumns, len(linearray)):
                                if delimiter in linearray[i]:
                                    linearray[i] = '"' + linearray[i] + '"'
                                datastring += delimiter + linearray[i]
            if refmismatch or extraalleles:
                if extraalleles:
                    fillin = 'Too many observed alleles'
                if refmismatch:
                    fillin = 'Mismatched reference'
                datastring = ''
                for i in range(0,headercolumns):
                    if datastring:
                        datastring += delimiter
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += linearray[i]
                for i in range(0,30):
                    datastring += delimiter + fillin
                for i in range (headercolumns, len(linearray)):
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += delimiter + linearray[i]
            elif not founddata:
                datastring = ''
                for i in range(0,headercolumns):
                    if datastring:
                        datastring += delimiter
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += linearray[i]
                for i in range(0,29):
                    datastring += delimiter + '0'
                datastring += delimiter + str(float((1/65000)*(1/65000)))
                for i in range (headercolumns, len(linearray)):
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += delimiter + linearray[i]
        else: #this is what we do if the current line being annotated is the same base as the last line (meaning we can skip over finding the data for it again)
            repeatlines += 1
            datastring = datastring.strip('\n')
            datastringarray = datastring.split(delimiter)
            linearrayindex = 0
            for i in range (0,len(datastringarray)):
                if i < headercolumns or i >= headercolumns + 30:
                    datastringarray[i] = linearray[linearrayindex]
                    linearrayindex += 1
            datastring = ''
            for datum in datastringarray:
                if datastring:
                    datastring += delimiter
                if delimiter in datum:
                    datum = '"' + datum + '"'
                datastring += datum
        datastring += '\n'
        try:
            outputfile.write(datastring)
        except:
            quit('Error writing to output file.')
        lastchromosome = chromosome
        lastposition = position
        lastreference = reference
        lastobserved = observed
        try:
            line = inputfile.readline()
        except:
            quit('Error reading from input file.')
    print ('Annotated ' + str(newlines) + ' unique loci and ' + str(repeatlines) + ' duplicated lines.' )

def main():
    import time  #this module lets us determine how long the run took (it is used only once at the very start and once at the very end of the program)
    starttime = time.time() #mark the start time
    job = checkargs()
    file = job[0]
    jobtype = job[1]
    if jobtype == 'librarysplit':
        if not librarysplit(file):
            quit('Error creating checksum hash for library.')
    elif jobtype == 'annotate':
        annotate(file)
    print ('Job completed successfully in ' + str(round(time.time() - starttime, 1)) + ' seconds.')
    
main()
    