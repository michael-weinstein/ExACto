#!/usr/bin/env python3
'''
Exome aggregation consortium annotation tool, Version 0.8
This program is designed to use the ExAC VCF data from the broad institute (downloadable from the link below).
ftp://ftp.broadinstitute.org/pub/ExAC_release
You will first need to run the library generation function using the commandline argument -s filename.vcf
After that, you will be able to annotate your variants using the commandline argument -f filename.txt
This version has not yet been tested against an ANNOVAR or VAX output, so it will probably encounter errors there.
Copyright 2014, Michael Weinstein, Cohn Lab, UCLA.  This is free to use for noncommercial entities, but please
let me know that you are using if I don't already know.  mweinste@ucla.edu
***THIS IS A BETA/TEST VERSION.  PLEASE REPORT BACK ANY UNEXPECTED BEHAVIOR OR ERRORS.  I HAVE TESTED THIS, BUT THERE MAY BE ERRORS I HAVE NOT YET FOUND.***
This version of ExACto requires using ExAC VCF version 0.2 or later
'''

def checkargs():  #subroutine for validating commandline arguments
    import argparse #loads the required library for reading the commandline
    import os  #imports the library we will need to check if the file exists
    file = False
    job = False
    parser = argparse.ArgumentParser()
    parser.add_argument ("-f", "--file", help = "Specify the desired file to annotate for submission.")  #tells the parser to look for -f and stuff after it and call that the filename
    parser.add_argument ("-s", "--split", help = "Specify a VCF to split into a subVCF library for use as a reference.")
    parser.add_argument ("-i", "--integrity", help = "Recreate the integrity check hash.", action = 'store_true')
    args = parser.parse_args()  #puts the arguments into the args object
    if  args.integrity:
        if args.file or args.split:
            usage("Rebuilding the integrity file can only be done on its own.  Please rerun witha  single command.")
            quit()
        print("Recreating the integrity file.")
        if createintegrity('subvcfs'):
            quit('Done.')
        else:
            quit('Unable to create a new integrity file.  Please verify that the directory \'subvcfs\' is a valid subdirectory of the folder containing this program.')
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
    outputindex = -1   #initializing the index counter for the output list to 0 (because lists start at 0 for the first element)
    for item in brokenlist:
        if outputlist and re.match('^".*[^"]$', outputlist[outputindex]):  #first checks to see if the list has been populated at all, then checks if the item in the current cell in the output list starts with quotes, but does not end with them
            outputlist[outputindex] = outputlist[outputindex] + delimiter + item  #append the next item from the broken, separated by a delimiter (fixes the broken entries)
        if not outputlist or not re.match('^".*[^"]$', outputlist[outputindex]): # if we are either just starting to populate the list or are looking at a complete item (something without unclosed quotes, even if we JUST completed them)
            if outputlist: #prevents us from trying to strip quotes off of list index -1, which cannot exist and trying to access it makes python very cross with us
                outputlist[outputindex] = outputlist[outputindex].strip('"') #removes the quotes from the entry in the list (they only get in the way at this point)  we will have to return them before writing them to a file
            outputindex += 1  #if the item in the output list looks complete
            outputlist.append(item) #create a new cell in the list containing the next item from the broken list (we will check to see if it's broken on the next iteration)
    return outputlist  #returns the list to the program for use
    
def getkeycolumns(headers):  #subroutine to extract a tupple containing the column numbers for key data points (location, reference, and observed alleles)
    import re
    chromosomefound = False  #initialize this to false
    positionfound = False  #initialize this to false as well, if we don't change them during the course of this subroutine, they will come back false and trigger an error message in the program
    referencefound = False
    observedfound = False
    for i in range (0, len(headers)):  #iterate through the list of headers
        if not chromosomefound and (re.search('^\W*?chr\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?chromosome\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?chrom\W*?$', headers[i], re.IGNORECASE)):  #use a regex to test if the column name indicates that it is a useful value to us (in this case chromosome), if so, we assign that value to the appropriate variable.  If we never find one, it remains false and returns false in that position to the annotation subroutine (causes the program to stop, as we need these values)
            chromosomecolumn = i
            chromosomefound = True
        elif not positionfound and (re.search('^\W*?start\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?position\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?pos\W*?$', headers[i], re.IGNORECASE)):
            positioncolumn = i
            positionfound = True
        elif not referencefound and (re.search('^\W*?ref\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?referencebase\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?Reference\W*?$', headers[i], re.IGNORECASE)):
            reference = i
            referencefound = True
        elif not observedfound and (re.search('^\W*?obs\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?samplealleles\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?alt\W*?$', headers[i], re.IGNORECASE) or re.search('^\W*?Alternate\W*?$', headers[i], re.IGNORECASE)):
            observed = i
            observedfound = True
    if not chromosomefound and positionfound and referencefound and observedfound:
        return False #if any of these values were not set, we do not return a tupple, instead we just return a false value
    return (chromosomecolumn, positioncolumn, reference, observed)  #if all these values were set, we return the tupple

def createdegeneratehash(): #not being called at the moment, but may become necessary later
    hashtable = {'A':['A'],'T':['T'],'G':['G'],'C':['C'],'N':['A','T','G','C'],'B':['T','G','C'],'D':['A','T','G'],'H':['A','T','C'],'K':['T','G'],'M':['A','C'],'R':['A','G'],'S':['G','C'],'V':['A','G','C'],'W':['A','T'],'Y':['C','T']}
    return hashtable  #line above and this line create a hashtable for degenerate bases with each letter in the hash connecting to a list of potential bases

def createzerohash(populations): #this subroutine returns a hash full of some set values (although it could easily be modified to return a hash full of anything) that will have the same structure as the hash of values when a variant is known and in the ExAC vcf
    datapoints = ['AC','AN','Hom'] #a list of the data points we need to fill in (I will go over details of the data structure for these hashes in the getfrequencyhash subroutine)
    fillins = {'AC':0,'AN':1,'Hom':0} #creates a dictionary of what to fill each data point's population:value frequencies with
    frequencyhash = {} #initializes an empty hash
    for point in datapoints: #iterate through each of the defined datapoints
        for population in populations: #iterate through each of the defined populations at each datapoint
            try: #the purpose of this try/except series is to avoid errors by putting a hash of a hash in a location that didn't already have an empty hash initialized to it (you can do that in PERL, but not python)
                frequencyhash[point][population] = fillins[point]  #if the empty hash for the datapoint already exists (or it's not an empty hash because there are already key:value pairs in it), write in the appropriate fillin
            except KeyError: #if trying to do the above raises a key error because the hash hasn't been initialized for the datapoint
                frequencyhash[point]={}  #first we create the empty hash
                frequencyhash[point][population] = fillins[point] #then we add a value to it
    return frequencyhash #and return the hash to the higher-level subroutine

def getfrequencyhash(line, populations, position):  #this is the "pointy end" of the annotation subroutine, this takes in the line from the ExAC reference and extracts the frequency data into a well-organized hash table

#Data structure
#hashname         frequencyhash
#                /      |      \
#datapoint     AC       AN      Hom
#         - - |        |        |
#         |   |-Pop:## |-Pop:## |-Pop:##
#Population   |-Pop:## |-Pop:## |-Pop:##
#         |   |-Pop:## |-Pop:## |-Pop:##
#         |   |-Pop:## |-Pop:## |-Pop:##
#         |   |-Pop:## |-Pop:## |-Pop:##
#         - - |-Pop:## |-Pop:## |-Pop:##
#any individual datapoint can be recalled as hashname[datapoint][population]

    import re
    multiallele = False
    populationset = ''
    frequencyhash = {}
    for population in populations: #iterate through populations
        if populationset: #if there is already a value written in the populationset string
            populationset += '|' #add in a pipe character
        populationset += population #add the current population to the string
    regex = re.compile('(AC|AN|Hom)_(' + populationset + ')=(.+?);')  #this regex is the "tip of the pointy end."  This regex will look for anything matching the datapoints (AC, AN, or Hom) followed by an underscore, then any of the populations (we built that part of the search using the populationset string), then the equal sign, then any number of characters over zero before a semicolon (we can't restrict to numbers because there will be a comma in there for multi-allele variants).  This will then return all those values as a list of tupples where the first value of each tupple is the datapoint, the second is the population, and the third is the actual number (or set of numbers for multi-allele)
    dataset = re.findall(regex, line) #this takes the regex we compiled above (we did not actually run it until this line) and finda ALL instances of it in the data line
    dataset = [list(entry) for entry in dataset]  #turns a tupple into a list so that it can be modified properly later (this is not a big issue until we have a multi-allele variant, a tupple won't allow for splitting the values within a cell)
    for datum in dataset: #iterate through the dataset (that 2-dimensional list of tupples, now lists, returned by the regex)
        if ',' in datum[2]: #check and see if there is a comma in the third element of each 3-element list as we iterate through them (a comma would indicate a locus with multiple reported alternate alleles)
            multiallele = True #set the value multiallele to true (we don't do anything to alter the data yet, we are just seeing if we need to)
            break #end the loop, once we find a comma, we know it's a multiallele case and don't need to waste cycles to confirm it repeatedly
    if multiallele: #if the previous test for multiallele loci came back positive
        for i in range (0,len(dataset)):  #iterate through the dataset
            if ',' in dataset[i][2]: #if the element being inspected contains a comma
                dataset[i][2] = dataset[i][2].split(',') #split that element into its own sublist on the comma
    for datum in dataset: #now iterate through the dataset again
        if multiallele and datum[0] != 'AN':  #iterate through the dataset list using the first element of each 3-element set as the datapoint key for the frequencyhash, the second element as the population, and the third as the value to set it equal to.
            try: #the use of try statements here, like in other cases where we are building a hash table, is to see if the key has already been initialized, and to initialize it if not
                frequencyhash[datum[0]][datum[1]] = datum[2][position] #note that we have to handle a multi-allele case differently: because some of the values are actually lists of values, we need an additional variable to locate our specific number (position within the list)
            except KeyError:
                frequencyhash[datum[0]]={}
                frequencyhash[datum[0]][datum[1]] = datum[2][position]
        else:
            try:
                frequencyhash[datum[0]][datum[1]] = datum[2] #this part of the subroutine does the same thing as the previous block, but does not have the position listed, as there will be only one element
            except KeyError:
                frequencyhash[datum[0]]={}
                frequencyhash[datum[0]][datum[1]] = datum[2]
    return frequencyhash

def addvcfinfolines(outputfile, numberalleles, datapoints, populations):
    poplong = {'AFR':'African/African American','AMR':'American','EAS':'East Asian','FIN':'Finnish','NFE':'Non-Finnish European','SAS':'South Asian','OTH':'Other'} #creates a dictionary for the population abbreviations
    datalong = {'AF':'Allele frequency (Allele observed/chromosomes counted)','HOMOF':'Homozygote Frequency (Homozygotes observed/[chromosomes counted/2])'} #creates a dictionary for the datapoint abbreviations
    for allele in range(0, numberalleles):
        for datapoint in datapoints:
            for population in populations:
                outputstring = ('## A' + str(allele +1) + '_' + datapoint + '_' + population + '\t: Allele' + str(allele +1) + ' ' + datalong[datapoint] + ' in ' + poplong[population] + ' from ExACto\n')
                try:
                    outputfile.write(outputstring)
                except:
                    quit('Error writing to output file.')
    try:
        outputfile.write('## F_rarest_allele\t: Highest population frequency of the rarest allele observed at locus from ExACto\n ## Combo_max\t : Maximum expected frequency of the allele combination from ExACto\n')
    except:
        quit('Error writing to output file.')


def annotate(file):  #and now for our main event
    import os
    import re
    if not os.path.isdir('subvcfs'): #if the library isn't already made or can't be found
        usage('No subvcf library detected in this directory.  Please be sure the script is running from the same directory that contains the library.')
        quit()
    fileformat = 'tdt'  #default format is tab-delimited
    delimiter = '\t'  #meaning that the delimiter is a tab
    lastlibraryfile = False  #initializes the lastlibrary file to false (we will use this to remember which library file we used in the last iteration so we don't reopen it if we don't have to)
    populations = ['AFR','AMR','EAS','FIN','NFE','SAS','OTH'] #initializing a list of the population in the ExAC
    datapoints = ['AF','HOMOF'] #initializing a list of the datapoins we will output for each population
    summarycolumns = ['F_rarest_allele','Combo_max']
    newcolumns = (((len(populations) +1) * len(datapoints)) * 2) + len(summarycolumns)  #logic here is that the number of columns we add is one for each datapoint for each population plus one for the maxvalue times two possible alleles for the line plus two for the columns we add that summarize population data (such as the frequency of the rarest allele and combo max) 
    if not checkintegrity('subvcfs'):  #checks the subvcf directory to be sure that no directories have been gained or lost or had their name changed since it was created
        quit('Check of subvcf library directory not passed.  Please be sure no files have been added, removed, or renamed in the library director')
    try:
        inputfile = open(file, 'r') #opens the file to be annotated
    except:
        quit('Error opening input file')
    outputfilename = file + '.ExACtoed.txt'  #creates the name of the output file
    if os.path.isfile(outputfilename):  #halts the program if there is already an annotated version of the file present  #off for debugging only  #this will end the program to avoid overwriting an existing annotation.  The user can delete or move it themselves if they want to rerun the annotation.
        usage('Output file already exists, please move, delete, or rename the existing output file for this source.')    #off for debugging only
        quit()                                                                                                           #off for debugging only
    try:
        outputfile = open(outputfilename, 'w', 1) #opens the output file
    except:
        quit('Error opening the output file.')
    try:  #This try/except (like most of the others you will see here) is to either read or write to a file, and end more gracefully if something fails in the process
        header = inputfile.readline().strip('\r\n') #reads the first line from the file to annotate.  This should be the header with all the column names.  We need to analyze this to determine where the important columns are
    except:
        quit('Error reading header from input file.')
    doublehashcount = 0
    vcflike = False
    vcfinfolineswritten = False
    while (header[0] == '#' and header[1] == '#') or not ((',' in header or '\t' in header) and re.search('chr', header, re.IGNORECASE) and re.search('ref', header, re.IGNORECASE)):  #this was added to handle VAX and other outputs that may have several double hashed lines before the headers and may have a hash starting the line (we will strip that away in a bit).  This looks for a line that is not double-hashed, contains either a comma or a tab delimiter, and has "CHR" and "REF" somewhere in it (case insensitive).  It keeps loading new lines until it finds one that fits the description.
        doublehashcount += 1
        if doublehashcount == 10:
            vcflike = True
        try:
            outputfile.write(header + '\n')
        except:
            quit('Error writing to output file')
        try:
            header = inputfile.readline().strip('\r\n')
        except:
            quit('Error reading header from input file.')
        if vcflike and not (header[0] == '#' and header[1] == '#') and not vcfinfolineswritten:
            addvcfinfolines(outputfile, 2, datapoints, populations)
            vcfinfolineswritten = True
    if not '\t' in header: #if we analyze the header line and find that there isn't a tab in there (indicating that tabs are not the delimiter)
        if not yesanswer('This does not appear to be a tab-delimited file.  Is it comma-delimited?'): #ask the user if this is a comma separated file
            quit('Please talk to Mike about getting this script to take your file format.')
        else:
            fileformat = 'csv'  #if the user says yes, the we treat this input as a CSV
    if fileformat == 'tdt': #probably unnecessary, as we already initialized the value, but good to be explicit here
        headerlist = header.split('\t')  #splits the header on each tab
        delimiter = '\t'
    if fileformat == 'csv':
        headerlist = quotedsplit(header) #if the file is a CSV, we split the header on commas.  This requires a full subroutine, since commas can be used in a field as well as being used to separate them and quoted fields make a huge difference
        delimiter = ','
    headercolumns = len(headerlist) #the variable headercolumns is set to the number of elements in the headerlist (this value is quite useful later)
    keycolumns = getkeycolumns(headerlist) #this subroutine takes the headerlist and returns the numbers of the columns containing values we will need in a tupple in a specific order
    if not keycolumns:  #if we could not get all those column numbers, we quit
        quit('Unable to find columns containing chromosome and/or position data.')
    chromosomecolumn = keycolumns[0] #this takes the values from the keycolumns tupple and gives them easier-to-remember names
    positioncolumn = keycolumns[1]
    referencecolumn = keycolumns[2]
    observedcolumn = keycolumns[3]
    try:
        line = inputfile.readline().strip('\r\n') #reads a new line from the user's annotation file (this should be the first line of data)
    except:
        quit('Error reading first data line from input file.')
    linenumber = 0  #initializes a counter of lines
    outputgroups = [] 
    for population in populations:
        outputgroups.append(population) #iterates through the populations and adds each one to the outputgroups list
    outputgroups.append('max') #then adds the max value to the list (we will calculate the max for each later)
    outputheaderlist = []
    for header in headerlist: #iterate through the headers that came with the file we are annotating
        outputheaderlist.append(header) #adds each header column name to the list of heaers to output
    for allele in range(0,2): #iterate through possible output alleles (in this case, just called A1 and A2)
        for point in datapoints: #iterate through the possible datapoints for each allele
            for group in outputgroups: #and then iterate through the different populations with 'max' as the last value (we added that in a few lines ago)
                outputheaderlist.append('A' + str(int(allele) + 1) + '_' + point + '_' + group) #add each of these values to the list of headers for the output file
    for summarycolumn in summarycolumns:
        outputheaderlist.append(summarycolumn)  #add each summary column (such as F_rarest_allele: frequency of rarest allele for the population in which it is most common)
    outputheaderstring = ''
    for item in outputheaderlist:  #iterate through each item in the outputheader list
        if outputheaderstring: #if we already have written values to the string of output headers (what we'll actually use to write to the file)
            outputheaderstring += delimiter #we add a delimiter to go between the values
        outputheaderstring += item #we add the item to the string
    outputheaderstring += '\n' #and cap the string with an end of line character
    try:
        outputfile.write(outputheaderstring) #then write the whole thing to the output file (which will now have its header line)
    except:
        quit('Error writing header to output file.')
    lastchromosome = False #initialize a bunch of values before starting a loop that will iterate through the datalines on the file
    lastposition = False
    lastreference = False
    lastobserved = False
    newlines = 0
    repeatlines = 0
    #reflineout = open('rfo.txt', 'w', 1)  #debugging only
    while line: #iterates through the file (so long as we load the next line to analyze before starting it).  If we have a blank line with nothing on it, this loop will end
        line = line.strip('\r\n') #take the new line and strip off any unwanted break or delimiter characters from the ends 
        linenumber += 1 #increment the line number counter
        #if linenumber >= 248: #debug only
        #    pass #debug only
        if line[0] == '#':  #if the first character of the line is a hashtag (indicating a comment and not data), 
            print('Skipped line ' + str(linenumber) + ' as it appears to be a comment.')  #tell the user we are skipping the line
            try:
                outputfile.write(line)  #write the line straight to the output file
            except:
                quit('Error writing to output file.')
            line = inputfile.readline() #read a new line
            continue #and restart the loop
        print ('Processing line ' + str(linenumber), end = '\r') #update the progress counter displayed to the user
        if fileformat == 'tdt': #this pair of IF statements just uses the filetype that we determined earlier to figure out how to split the lines (as above, if the file is CSV, it may have quoted elements and those need to be handled differently)
            linearray = line.split('\t')
        elif fileformat == 'csv':
            linearray = quotedsplit(line)
        chromosome = re.sub('^chr', '',linearray[chromosomecolumn], flags=re.IGNORECASE) #because we want only the name of the chromosome, but sometimes chromosomes are recorded as chr##, this will remove chr from the beginning of the chromosome value (and do so regardless of the case of the letters).  Then it will save the value as chromosome
        position = linearray[positioncolumn] #this saves the position value from the data line to position (this is simple and we don't have to worry that someone put anything in front of it).  The next two lines do the same thing
        reference = linearray[referencecolumn]
        observed = linearray[observedcolumn]
        if not (chromosome == lastchromosome and position == lastposition and reference == lastreference and observed == lastobserved): #This asks if the current chromosome, position, reference, and observed alleles are the same as the last time we iterated through.  In an output listing several possible transcripts for each variant, this recycling of already-calculated values cuts more than 50% of the time required.
            newlines += 1  #if not, we increment the count of newlines (this is only for display to the user at the end and not used for analysis)
            frequencyhash = {}  #initializing a bunch of values
            refmismatch = False
            refandaltswapped = False
            extraalleles = False
            observedarray = []
            observed = observed.strip('"')
            if '/' in observed: #series of if statements just checks for any splitting characters (comma or slash) in the observed alleles field.  If so, it splits it into a list on the character.  If not, it just creates a one element table from it (so that we always have a list, and don't have to worry about ending up with a string for a single allele)
                observedarray = observed.split('/')
            elif ',' in observed:
                observedarray = observed.split(',')
            else:
                observedarray.append(observed)
            currentlibraryfile = str(chromosome) + 'c' + str(int(position) - (int(position) % 100000)) + '.subvcf'  #uses the chromosome and position data to name the library file to look for ExAC reference data
            if os.path.isfile('subvcfs/' + currentlibraryfile): #checks if the library file exists
                founddata = True #if so, sets this value appropriately
            else:
                founddata = False  #if not, sets it to false so we can stop looking
            if founddata: #if we did find the appropriate library file
                if currentlibraryfile != lastlibraryfile:  #if it is not the same library file as last time, we open a new one, read the whole thing into memory (making regex searches faster) and then close it
                    try:
                        libraryfile = open('subvcfs/' + currentlibraryfile, 'r')
                        library = libraryfile.read() #slurps the whole file into the library string
                        libraryfile.close()
                    except:
                        quit('Error opening subvcf library file.')
                refline = re.search('^(' + str(chromosome) + '\t' + str(position) + '\t.*?)$', library, re.MULTILINE) #search the library for a line with the same chromosome and position, then take that entire line and save it as refline
                if not refline: #if refline is blank because we did not find anything
                    founddata = False  #set this value so the rest of this program knows to treat it as an unseen variant
                else:
                    refline = refline.group(0) #otherwise we extract the line from the regular expression returned value (and get rid of the rest of the regex junk)
                    #reflineout.write(refline + '\n')  #debugging only
                    reflinearray = refline.split('\t')  #and split the refline on each tab (we know that's its delimiter)
                if founddata: #if we found data for that locus (if we didn't, we get to skip all this work)
                    exacalts = reflinearray[4].split(',') #creates an array of alternate alleles seen in exac
                    exacalthash = {}
                    for i in range (0, len(exacalts)): #creates a hash with the allele number:allele sequence as pairs
                        exacalthash[exacalts[i]]= i
                    if refline and (reflinearray[3] != reference):  #makes sure again that we have some value in refline (we shouldn't be here if we don't) and then checks to make sure that the variant caller called the same reference allele as ExAC has
                        refmismatch = True  #set this value so the rest of the program knows the line has mismatched references
                        if reflinearray[3] in observedarray:  #if the reference value from ExAC was observed as an alternate allele in the data
                            for allele in observedarray:  #and if we iterate through each allele
                                if allele == reflinearray[3] or allele in exacalts:  #and find that each allele was listed as either the reference or alternate allele in ExAC
                                    refandaltswapped = True  #we set this value to true, since the reference and alternates were swapped between sources.  We will still not make a call, but the message displayed will be different (this sometimes happens with very common variants where differnt sources think one is more common and call it reference)
                        if refandaltswapped:
                            print('Reference and alternate alleles swapped between input and ExAC on line ' + str(linenumber) + '.')
                        else:
                            print('Reference base mismatch starting on line ' + str(linenumber) + '.')  #let the user know there was a mismatch (this is important because a bunch of mismatches can indicate different Hg versions used to annotate between the reference and variant files)    
                    else:
                        try:
                            observedarray.remove(reference)  #if the observed alleles contains the reference allele, we remove it (since we don't analyze frequency data from it).  At some point, we could probably deduce the data for it from the listed data, and then we could include that
                            if len(observedarray) > 1: #if we do this and still have more than one allele left, we have three reported alleles, and need to report the line as an error (possibly a bad read or mosaic)
                                extraalleles = True
                            singlenonref = True  #if we could remove that value without error, then we have a single non-reference allele
                        except:
                            singlenonref = False  #if not, then both alleles are non-reference
                            if len(observedarray) == 1: #if only one observed allele is listed, it is homozygous and non-reference
                                homozygousrare = True
                            elif len(observedarray) > 1 and observedarray[0] == observedarray[1]: #if we have more than one element in the observed array and those elements are the same
                                homozygousrare = True #the subject is homozygous for a rare allele
                                observedarray.remove(observedarray[0]) #and we remove the duplicated allele so we don't analyze it twice
                                if len(observedarray) > 1: #if we remove one allele from an identical pair and still have two alleles on the set, there was a third allele reported
                                    extraalleles = True
                            else:
                                if not extraalleles: #if we did not get back an indication of three alleles or more at the locus
                                    homozygousrare = False #the subject is heterozygous at the locus for two non-reference alleles
                        if not extraalleles and len(observedarray) in range(1,3):
                            extraalleles = False #one last check to make sure we don't have extra alleles
                            for allele in (observedarray): #go through each of the observed alleles
                                try: #use these try/except statements to build a frequency hash for each allele (the frequency hash data structure is better explained in the getfrequencyhash subroutine)
                                    frequencyhash[allele] = getfrequencyhash(reflinearray[7], populations, exacalthash[allele])  #the frequency hash for each allele is set to the returned value (a hash of hashes itself) from the getfrequencyhash subroutine
                                except KeyError:  #if it returns a keyerror (because the hash of alternate alleles observed for ExAC did not list this allele), then this variant is unique and we need a hash full of zeroes for the frequency (and ones for the chromosome count, to avoid divisio by zero)
                                    frequencyhash[allele] = createzerohash(populations)
                        else:
                            print('Possible error on line ' + str(linenumber) + ', subject appears to have more than 2 alleles.') #if we hit this statement because there were extra alleles on the line (more than 2)
                            extraalleles = True #set this value (which should already be set)
                        if not extraalleles: #otherwise, we continue the analysis
                            i = 0  #initialize a bunch of values
                            namelist = []
                            namehash = {}
                            valuearray = []
                            for allele in observedarray: #iterate through the alleles in the observed array
                                for point in datapoints: #iterate throught the datapoints
                                    for population in populations: #iterate through the populations
                                        try: #use this try/except stack to build up a hash containing all of our frequency values. The value for each one will be its position in a list containing only location values
                                            namehash[allele][point][population] = i
                                        except KeyError:
                                            try:  #we do this if we get a key error indicating that we don't yet have a hash declared for the datapoint
                                                namehash[allele][point] = {}
                                                namehash[allele][point][population] = i
                                            except KeyError: #and we do this if we still get a key error, indicating that we don't have a hash declared for the allele itself
                                                try:
                                                    namehash[allele] = {}
                                                    namehash[allele][point] = {}
                                                    namehash[allele][point][population] = i
                                                except:
                                                    print('This is what you get for not making the data structure simpler')
                                        namelist.append([allele, point, population])  #add to the list of field names
                                        if point == 'AF':  #if we are analyzing allele frequency, that equates to the number of times seen in the population divided by the chromosomes counted for the population
                                            try:
                                                entry = int(frequencyhash[allele]['AC'][population])/int(frequencyhash[allele]['AN'][population])
                                            except KeyError:  #if we get a key error, that means that the allele is unique, so the frequency for it being seen is zero
                                                entry = 0
                                            except ZeroDivisionError:  #if we get a zero division error, that means the locus had no coverage in the population and we cannot make a call
                                                entry = 'NA'
                                        elif point == 'HOMOF':  #if we are analyzing homozygous frequency for the allele
                                            if chromosome not in ('X','x','Y','y','MT','Mt','mt'): #first we check and make sure it's not a sex chromosome (due to issues with hemizogosity) or the mitochondrial DNA
                                                try:
                                                    entry = int(frequencyhash[allele]['Hom'][population])/(int(frequencyhash[allele]['AN'][population])/2)  #if it is an autosome, homozygosity frequency is estimated as homozygotes counted/(chromosomes/2)
                                                except KeyError:  #if we get a key error (meaning the allele is not reported in exac), we assume the allele to be unique (or at least have zero frequency in this database) and have zero homozygotes
                                                    entry = 0 
                                                except ZeroDivisionError: #and if we get no coverage at the locus, we return NA
                                                    entry = 'NA'
                                            else:
                                                entry = 'NA'  #we also return NA for non-autosomal variants here to avoid giving false information in the case of hemizygosity (GitHub currently reports that ExAC does not handle hemizygosity for sex chromosomes very well)
                                        valuearray.append(entry)
                                        i += 1
                                    maximum = 0
                                    namehash[allele][point]['max'] = i  #create an entry in the name hash for max (the maximum frequency in the populations)
                                    namelist.append([allele, point, 'max']) #add this to the name list as well
                                    for population in populations: #iterate throught the populations
                                        if valuearray[namehash[allele][point][population]] != 'NA':  #if the value was not reported as NA (indicating no coverage)
                                            if maximum < valuearray[namehash[allele][point][population]]: #if the current value being inspected is greater than the previous maximum
                                                maximum = valuearray[namehash[allele][point][population]] #that value becomes the new maximum
                                    valuearray.append(maximum) #after iterating through all populations, we have the maximum and that gets entered in the array of values
                                    i += 1 #and we increment the index pointer to the next cell in the list
                            if not singlenonref and not homozygousrare:  #if this allele is heterozygous for two non-reference alleles
                                allele1max = valuearray[namehash[observedarray[0]][point]['max']] #this just looks up the maximum frequency for allele 1
                                allele2max = valuearray[namehash[observedarray[1]][point]['max']] #this does the same for allele 2
                                rarestallele = min(int(allele1max),int(allele2max)) #this figures out which one is less
                                if allele1max == 0:
                                    allele1max = 1/65000    #this changes a zero value for these into 1/65000 to give us a more reasonable estimate of combination frequency (note that this will still treat a unique allele as being exceedingly rare)
                                if allele2max == 0:
                                    allele2max = 1/65000
                                combomax = allele1max * allele2max
                            elif singlenonref:
                                rarestallele = valuearray[namehash[observedarray[0]]['AF']['max']]  #if only 1 allele was non-reference, we just use its value for the rarest and don't worry about the combo
                                combomax = 'NA'
                            elif homozygousrare:
                                rarestallele = valuearray[namehash[observedarray[0]]['AF']['max']]  #and if the person is homozygous for a non-reference allele, the expected frequency of the combination (a homozygote) would be the highest frequency squared
                                if rarestallele == 0:
                                    combomax = (1/65000) ** 2
                                else:
                                    combomax = rarestallele ** 2
                            datastring = ''  #initialize an empty string for building our output line
                            for i in range(0,headercolumns):  #iterate through headercolumns (we use an i and a range so we can have a pointer of our position).  We are only going to the end of our named columns from the header (some outputs have additional data at the end of the line in un-headed columns; we will get to those)
                                if datastring:  #if there is already something written to the string, we add a delimiter before adding anything else
                                    datastring += delimiter 
                                if delimiter in linearray[i]:  #if the currently-used delimiter character is in the data we want to write (usually a comma in a CSV)
                                    linearray[i] = '"' + linearray[i] + '"'  #we add quotes to the beginning and end of the value to make sure that it is kept together
                                datastring += linearray[i]  #then add the value to the end of the string
                            for value in valuearray: #now iterate through the new values we want to add and put those on the string (the string is already started, so we don't have to worry about not sticking on a delimiter before the first item)
                                datastring += delimiter + str(value)
                            if len(observedarray) == 1: #this is a little confusing, but if we only had a single non-reference allele, regardless of zygosity, we only have one set of frequencies to report.  This fills in the Allele 2 columns with NA
                                for value in valuearray:
                                    datastring += delimiter + 'NA'
                                #if homozygousrare:  #if the allele is homozygous for a nonreference
                                #    rarestallele = valuearray[namehash[observedarray[0]][point]['max']] #this just looks up the maximum frequency for allele 1
                                #    combomax = rarestallele ** 2
                                #elif singlenonref: #if the locus is heterozygous for reference and nonreference
                                #    rarestallele = valuearray[namehash[observedarray[0]][point]['max']] #this just looks up the maximum frequency for allele 1
                                #    combomax = 'NA'
                            datastring += delimiter + str(rarestallele) + delimiter + str(combomax)  #add the values to the datastring for output
                            for i in range (headercolumns, len(linearray)):  #iterate through the columns from the original annotation that came after the ones with headers (in other words, headerless columns at the end of the table)
                                if delimiter in linearray[i]: #as before, if the element we want to add has our delimiter in it, we need to quote the element to keep it together
                                    linearray[i] = '"' + linearray[i] + '"' 
                                datastring += delimiter + linearray[i] #and then add it to the growing string for eventual output
            if refmismatch or extraalleles: #This handles what to output if there was a problem with the line (reference mismatch or 3+ alleles).  We will fill in the values with a message indicating why we did not give a value
                if extraalleles:
                    fillin = 'Too many observed alleles'
                if refmismatch: #if reference was mismatched between the exome data and the ExAC data, we want to tell if the values were swapped (indicating a likely common allele where the reference has changed between genome edits) or if the two were entirely different
                    if refandaltswapped:
                        fillin = 'Swapped reference and alternate alleles (likely common).'
                    else:
                        fillin = 'Mismatched reference'
                datastring = ''
                for i in range(0,headercolumns):
                    if datastring:
                        datastring += delimiter
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += linearray[i]
                for i in range(0,newcolumns):  #changed for number of columns added
                    datastring += delimiter + fillin
                for i in range (headercolumns, len(linearray)):
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += delimiter + linearray[i]
            elif not founddata:  #if we found no ExAC reference data, we treat the variant as unique and put out the appropriate values
                datastring = ''
                for i in range(0,headercolumns):
                    if datastring:
                        datastring += delimiter
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += linearray[i]
                for i in range(0,int((newcolumns-len(summarycolumns))/2)):  #we have to retype as int here because division forces the value to a float type, even though it should always be the division of an even number by 2
                    datastring += delimiter + '0'
                for i in range(0,int((newcolumns-len(summarycolumns))/2)):
                    datastring += delimiter + 'NA'
                datastring += delimiter + '0'
                datastring += delimiter + str(float((1/65000)*(1/65000)))
                for i in range (headercolumns, len(linearray)):
                    if delimiter in linearray[i]:
                        linearray[i] = '"' + linearray[i] + '"'
                    datastring += delimiter + linearray[i]
        else: #this is what we do if the current line being annotated is the same base as the last line (meaning we can skip over finding the data for it again)
            repeatlines += 1 #add one to the number of repeated lines we ran
            datastring = datastring.strip('\n')  #take the previous iteration's datastring and remove the end of line on the end
            datastringarray = datastring.split(delimiter) #split the datastring back into a table on the delimiter
            linearrayindex = 0  #this keeps track of where in the linearray we are looking
            for i in range (0,len(datastringarray)): #iterate through the datastring array
                if i < headercolumns or i >= headercolumns + newcolumns: #indexing here depends on how many columns we add, double check if problems #if we are looking at a column that was not generated by this program (since values generated here will not change for different transcripts of the same gene, and it saves time not to calculate them again), change each value in the data string array to reflect the current line (most of these will actually remain the same, except for things pertaining to the specific transcript)
                    datastringarray[i] = linearray[linearrayindex]
                    linearrayindex += 1 #increment the counter so we know we are
            datastring = '' #reinitialize the datastring (which is still going to contain the last iteration's data)
            for datum in datastringarray:  #then iterate through our datastringarray with some new values from the currentline and form a new datastring to output
                if datastring:
                    datastring += delimiter
                if delimiter in datum:
                    datum = '"' + datum + '"'
                datastring += datum
        datastring += '\n' #and cap off the string with an end of line
        try:
            outputfile.write(datastring) #and, regardless of which statement was used to write the datastring (new line or repeat line), we write it to the output file
        except:
            quit('Error writing to output file.')
        lastchromosome = chromosome  #these values remember this iteration's important information so that we can check next time to see if this is a repeat line or a new one
        lastposition = position
        lastreference = reference
        lastobserved = observed
        try:
            line = inputfile.readline()  #then we read a new line from the input file (if it reads blank becaue it's the end of the file, the loop will exit)
        except:
            quit('Error reading from input file.')
    print ('Annotated ' + str(newlines) + ' unique loci and ' + str(repeatlines) + ' duplicated lines.' )  #Tells the user a summary of what was done.
    #reflineout.close()  #debugging only

def main():
    import time  #this module lets us determine how long the run took (it is used only once at the very start and once at the very end of the program)
    starttime = time.time() #mark the start time
    job = checkargs() #make sure the commandline arguments were valid, if they were not, this subroutine will exit the program
    file = job[0] 
    jobtype = job[1]
    if jobtype == 'librarysplit':  #this if/elseif statement will tell the program to either generate a library or annotate a file passed in the argument
        if not librarysplit(file):
            quit('Error creating checksum hash for library.')
    elif jobtype == 'annotate':
        annotate(file)
    print ('Job completed successfully in ' + str(round(time.time() - starttime, 1)) + ' seconds.') #report back to the user how long this all took
    
main()
    