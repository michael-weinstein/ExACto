#!/usr/bin/env python3

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

def directoryhash(directory):
    import os
    import hashlib
    allfiles = os.listdir(directory)
    filestring = ''.join(allfiles)
    filestringhash = hashlib.md5(rawbytes(filestring))
    return filestringhash

def rawbytes(stringin):  #subroutine to take a string and make it into UTF-8 bytes, often for sending via socket connection
    return bytes(stringin, 'utf-8')

def checkintegrity(directory):
    hashfile = open('hashsum', r)
    recordedhash = hashfile.readline()
    hashfile.close()
    if recordedhash == directoryhash(directory):
        return True
    else:
        return False
        
def createintegrity(directory):
    hashfile = open('hashsum', 'w')  #make sure this is actually created and not just buffered at this stage.  close and reopen if needed
    hashfile.write(directoryhash(directory))
    hashfile.close()
    return True

def librarysplit(filename):
    import re
    import os
    if os.path.exists('subvcfs'):
        usage('Library directory already exists.  Please remove the old one before creating a new one')
        quit()
    try:
        os.makedirs('subvcfs')  #actually creates the new directory
    except PermissionError:
        print ('\nPermission to create the library directory was denied.  Please run this program with appropriate permissions or create it using your file manager.')
        quit()
    except:
        print ('\nUnable to make the library directory.')
        quit()    
    lastchromosome = 0
    linecount = 0
    libraryfile = False
    lastpositionblock = 0
    exac = open(filename, 'r')
    line = exac.readline()
    while line:
        linecount += 1
        if re.match('^\#', line):
            line = exac.readline()
            continue
        linearray = line.split('\t')
        chromosome = linearray[0]
        position = linearray[1]
        positionblock = str(int(position) - (int(position) % 1000000))
        if positionblock != lastpositionblock or chromosome != lastchromosome:
            if libraryfile:
                libraryfile.close()
            libraryfile = createlibraryfile(chromosome, positionblock)
        libraryfile.write(line)
        line = exac.readline()
        lastpositionblock = positionblock
        lastchromosome = chromosome
        print('Processed ' + str(linecount) + ' lines.', end = '\r')
    print('\n')
    createintegrity('subvcfs')
    return True

def createlibraryfile(chromosome, positionblock):
    libraryfilename = str('subvcfs/' + chromosome + 'c' + positionblock + '.subvcf')
    libraryfile = open(libraryfilename, 'w') #needs a try/except block
    return libraryfile

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

def quotedsplit(input, delimiter = ','):
    import re
    brokenlist = input.split(delimiter)
    outputlist = []
    outputindex = -1
    for item in brokenlist:
        if re.match('^".*[^"]$', outputlist[outputindex]):
            outputlist[outputindex] = outputlist[outputindex] + ',' + item
        else:
            outputindex += 1
            outputlist.append(item)
    return outputlist
    
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

def getfrequencyhash(line, populations, position):
    multiallele = False
    populationset = ''
    frequencyhash = {}
    for population in populations:
        if populationset:
            populationset += '|'
        populationset += population
    regex = re.compile('(AC|AN|Hom)_(' + populationset + ')=(.+?);')
    dataset = re.findall(regex, line)
    dataset = [list(entry) for entry in dataset]
    for datum in dataset:
        if ',' in datum[2]:
            multiallele = True
            break
    if multiallele:
        for i in range (0,len(dataset)):
            if ',' in dataset[i][2]:
                dataset[i][2] = dataset[i][2].split(',')
            else:
                dataset[i][2] = [dataset[i][2],dataset[i][2]]
    for datum in dataset:
        try:
            frequencyhash[datum[0]][datum[1]] = datum[2][position]
        except KeyError:
            frequencyhash[datum[0]]={}
            frequencyhash[datum[0]][datum[1]] = datum[2][position]
    return frequencyhash

def annotate(file):
    import os
    import re
    fileformat = 'tdt'
    delimiter = '\t'
    lastlibraryfile = False
    populations = ['AFR','AMR','EAS','FIN','NFE','SAS']
    datapoints = ['AF','HOMOF']
    if not checkintegrity(subvcfs):
        quit('Check of subvcf library directory not passed.  Please be sure no files have been added, removed, or renamed in the library director')
    inputfile = open(file, 'r')
    header = inputfile.readline()
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
    if not locuscolumns:
        quit('Unable to find columns containing chromosome and/or position data.')
    chromosomecolumn = locuscolumns[0]
    positioncolumn = locuscolumns[1]
    referencecolumn = locuscolumns[2]
    observedcolumn = locuscolumns[3]
    outputfilename = file + '.ExACtoed.txt'
    outputfile = open(outputfilename, 'w', 1)
    line = inputfile.readline()
    linenumber = 0
    outputgroups = populations
    outputgroups.append('max')
    outputheaderlist = headerlist
    for allele in range(0,2):
        for point in datapoints:
            for group in outputgroups:
                outputheaderlist.append('A' + int(allele + 1) + '_' + point + '_' + group)
    outputheaderlist.append('F_rarest_allele')
    outputheaderlist.append('Combo_max')
    outputheaderstring = ''
    for item in outputheaderlist:
        if outputheaderstring:
            outputheaderstring += delimiter
        outputheaderstring += item
    outputheaderstring += '\n'
    outputfile.write(outputheaderstring)
    while line:
        refmismatch = False
        lineerror = False
        linenumber += 1
        print ('Processing line ' + str(linenumber), end = '\r')
        if fileformat == 'tdt':
            linearray = line.split('\t')
        if fileformat == 'csv':
            linearray = quotedsplit(line)
        chromosome = linearray[chromosomecolumn]
        position = linearray[positioncolumn]
        reference = linearray[referencecolumn]
        observed = linearray[observedcolumn]
        observedarray = observed.split('/')
        currentlibraryfile = str(chromosome) + 'c' +str(int(position) - int((position % 1000000))) + '.subvcf'
        if currentlibraryfile != lastlibraryfile:
            libraryfile = open('subvcfs/' + currentlibraryfile, 'r')
            library = libraryfile.read() #slurps the whole file into the library string
            libraryfile.close()
        regex = re.compile('\n' + str(chromosome) + '\t' + str(position) + '\t.*\n')
        refline = re.search(regex, library, re.IGNORECASE)
        if not refline:
            frequencyhash = False
        else:
            reflinearray = refline.split('\t')
        if refline and (reflinearray[3] != reference):  #if this doesn't analyze lazy, we need two if statements
            print('Reference base mismatch on line ' + str(linenumber) + '.')
            refmismatch = True
        else:
            try:
                observedarray.remove(reference)
                singlenonref = True
            except:
                singlenonref = False
                if observedarray[0] == observedarray[1]:
                    homozygousrare = True
                    observedarray.remove(observedarray[0])
                else:
                    homozygousrare = False
            if len(observedarray) in range(1,3):
                for i in range(0,2):
                    frequencyhash[observedarray[i]] = getfrequencyhash(reflinearray[7], populations, i)
            else:
                print('Possible error on line ' + str(linenumber) + ', subject appears to have more than 2 alleles.')
                lineerror = True
        i = 0
        namelist = []
        namehash = {}
        valuearray = []
        for allele in observedarray:
            for value in ('AF','HOMOF'):
                for population in populations:
                    try:
                        namehash[allele][value][population] = i
                    except KeyError:
                        try:
                            namehash[allele][value] = {}
                            namehash[allele][value][population] = i
                        except KeyError:
                            try:
                                namehash[allele] = {}
                                namehash[allele][value] = {}
                                namehash[allele][value][population] = i
                            except:
                                print('This is what you get for not making the data structure simpler')
                    namelist.append([allele, value, population])
                    if value == 'AF':
                        entry = int(frequencyhash[allele]['AC_'+population])/int(frequencyhash[allele]['AN_'+population])
                    elif value == 'HOMOF':
                        if chromosome not in ('X','x','Y','y','MT','Mt','mt'):
                            entry = int(frequencyhash[allele]['Hom_'+population])/(int(frequencyhash[allele]['AN_'+population])/2)
                        else:
                            entry = 'NA'
                    valuearray.append(entry)
                    i += 1
                maximum = 0
                namehash[allele][value]['max'] = i
                namelist.append([allele, value, 'max'])
                for population in populations:
                    if maximum > valuearray[namehash][allele][population]:
                        maximum = valuearray[namehash][allele][population]
                valuearray.append(maximum)
                i += 1
        if not singlenonref and not homozygousrare:
            allele1max = valuearray[namehash][observedarray[0]]['AF']['max']
            allele2max = valuearray[namehash][observedarray[1]]['AF']['max']
            rarestallele = min(int(allele1max),int(allele2max))
            if allele1max == 0:
                allele1max = 1/65000
            if allele2max == 0:
                allele2max = 1/65000
            combomax = allele1max * allele2max
        elif singlenonref:
            rarestallele = valuearray[namehash][observedarray[0]]['AF']['max']
            combomax = 'NA'
        elif homozygousrare:
            rarestallele = valuearray[namehash][observedarray[0]]['AF']['max']
            combomax = (1/65000)^2
    datastring = ''
    for i in range(0,headercolumns):
        if datastring:
            datastring += delimiter
        datastring += linearray[i]
    for allele in range(0,2):
        for point in datapoints:
            for group in outputgroups:
                datastring += delimiter + frequencyhash[allele][point][group]
    for i in range (headercolumns, len(linearray)):
        datastring += delimiter + linearray[i]
    outputfile.write(datastring)
    line = inputfile.readline()

#create line to write for either output format
#???
#profit
                    
            
        
        
            
def main():
    import time  #this module lets us determine how long the run took (it is used only once at the very start and once at the very end of the program)
    starttime = time.time() #mark the start time
    job = checkargs()
    file = job[0]
    jobtype = job[1]
    if jobtype == 'librarysplit':
        librarysplit(file)
    elif jobtype == 'annotate':
        annotate(file)
    print ('\nJob completed successfully in ' + str(round(time.time() - starttime, 1)) + ' seconds.')
    
main()
    