
# coding: utf-8

# In[10]:


#OrfScanner

#importing libraries
import re
import argparse

#creating command line options
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="Identify the output file", type=str, default = 'orfs.fasta')
parser.add_argument("-r", "--readingframe", help="Limit the ORF search to 1 frame: 1, 2, 3, -1, -2 or -3", type=int, default = 0)
parser.add_argument("-c", "--cutoff", help="Set the minimum length of ORFs outputted - amino acid length", type=int, default = 100)
parser.add_argument("-f", "--inputfile", help="Input file, must be in FASTA format", type=str, nargs='+')
args = parser.parse_args()

#dealing with unexpected command line arguements
if args.readingframe > 3:
    print('Frame must be one of 1 2 3 -1 -2 -3. Default finds ORFs in all frames \n Exiting')
    exit()
if args.readingframe < -3:
    print('Frame must be one of 1 2 3 -1 -2 -3. Default finds ORFs in all frames \n Exiting')
    exit()
if args.cutoff <2:
    print('Minimum ORF length must be greater than 2 \n Exiting')
    exit()
if args.inputfile == None:
    print('Please enter at least 1 input file using the -f option \n Exiting')
    exit()

#fasta file input which can deal with more than 1 fasta file
filelist = list(args.inputfile)
#print("Specify at least 1 input file and then type 'run' to begin")
#a = str(input('Fasta file / run '))
#while a != 'run':
 #   filelist.append(a)
 #   a = str(input('Fasta file / run '))
#print(filelist)

#reading in each file
filecount = 0
for file in filelist:
    filecount += 1
    dna = str(file)
    cutoff =  args.cutoff
    if len(filelist) > 1:
        filename = str(filecount) + str(args.output)
    if len(filelist) == 1:
        filename = str(args.output)    
    fileout = open(filename,'w')
    fileObj = open(dna, 'r')
    linelist = fileObj.readlines() # linelist is a list of strings
    #finding the organisms name
    firstline = linelist[0]
    organism = firstline.replace('\n','') #removing line breaks
    #checking file is in fasta format based on > in first line and allows this to be overidden by the user
    if organism[0] != '>':
        print('Are you sure {0} is in FASTA format? No > on 1st file line'.format(file))
        ignore = input('Continue? Y/N: ')
        if ignore != 'Y':
            print('Exiting')
            exit()
    lines = linelist[1:]
    #mergeing rest of file into 1 string of nucleotides
    sequence_string = '' #merging lines in the list into a string
    for item in lines:
        sequence_string +=item

    sequence = sequence_string.replace('\n','') #removing line breaks

    #Creating the complementary strand using reverse complement
    def revComp (x): #creating function revComp
        """
        revComp creates the reverse complement of a nucelotide sequence.
        Parameters
        ----------
        arg1 : str
            A capitalized nucleotide sequence
        Returns
        -------
        str
            A reverse complement strand
        """
        transTab = str.maketrans('atgc', 'tacg')
        comp = x.translate(transTab)
        rev = comp[::-1]
        return rev 

    rc = revComp(sequence) #rc is the reversed complementary sequence
    
    #capitalizing sequences
    capsequence = sequence.upper()
    caprc = rc.upper()
    
    #translating the sequence and finding ORFs
    def orfsearch (aaseq = str, cutoff = int, readingframe = int):
        """
        orfsearch translates a nucleotide sequence to an amino acid sequence.
        It then searches for opening reading frames - non-nested sequences which begin with a start codon (M only) and end with a stop codon (*).
        orfsearch is able to search within all 3 reading frames of a strand.
        Parameters
        ----------
        aaseq : str
            A capitalized nucleotide sequence - only one strand, forward or reverse, to be supplied at once
        cutoff: int
            Minimum length of ORFs outputted - in amino acids
        readingframe: int
            1 , 2 or 3; determines which frame the nucelotide sequence is translated in
        Returns
        -------
        str
            Writes to a file all the ORFs found and gives the organism name, frame, ORF number and ORF position in the amino acid sequence followed by the amino acid sequence in the ORF.
            For example:
            >I.claudius Frame:1 ORF:1 Position:1
            ORFSEQUENCE
        """
        
        # user friendly naming of frames
        aminoseq = ''
        if readingframe == 0 and aaseq==caprc:
            framelabel = '-1'
        if readingframe==1 and aaseq==caprc:
            framelabel = '-2'
        if readingframe==2 and aaseq==caprc:
            framelabel = '-3'
        if readingframe==0 and aaseq==capsequence:
            framelabel = '1'
        if readingframe==1 and aaseq==capsequence:
            framelabel = '2'
        if readingframe==2 and aaseq==capsequence:
            framelabel = '3'
        
        #dictionary containing the genetic code - bacteria
        geneticcode = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT" : "N", "ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T", "AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S", "ATA" : "I", "ATC" : "I", "ATG" : "M", "ATT" : "I", "CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H", "CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P", "CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R",	 "CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L", "GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D", "GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A", "GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G", "GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V", "TAA" : "*", "TAC" : "Y", "TAG" : "*", "TAT" : "Y", "TCA" : "S", "TCC" : "S", "TCG" : "S", "TCT" : "S", "TGA" : "*", "TGC" : "C", "TGG" : "W", "TGT" : "C", "TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT" : "F" }
        
        #translating using the genetic code dicitonary
        for i in range(readingframe, len(aaseq), 3): #readingframe offsets the start of the iteration
            codon = str(aaseq[i:i+3]) 
            trans = geneticcode.get(codon, 'X') #if a gap or unknown base is encountered, X is put into the translated sequence
            aminoseq += trans
        
        #orf finding
        gflag=False #gflag will be used to prevent nested ORFs being taken - gflag = True when it's within an ORF
        g=''
        genelist =[]    
        position = 0
        frame = readingframe
        poslist = []
        
        #iterating over the amino acids in the translated sequence looking for ORFs and saving to a list
        for aa in aminoseq:
            position +=1
            if gflag==True and aa != '*':
                g+=aa
            if aa in 'M':                                                                          
                if gflag==False:
                    startpos = position
                    g+=aa
                gflag=True 
            if aa in '*':                                                                           
                gflag=False
                length = len(g)
                if length > cutoff:
                    genelist.append(g)
                    poslist.append(startpos)
                g =''
                count = -1

        #removing orfs which contained a gap or an unknown base - these were translated to X
        for gene in genelist:
            count += 1
            for amino in gene:
                if amino in 'X':
                    del(genelist[count])
                    del(poslist[count])
        #user message to show number of ORFs found
        print('{0}: {1} ORFs found in frame {2}'.format(str(file), str(len(poslist)), str(framelabel)))
        
        #writing output to the output file
        count2 = -1
        for i in genelist:
            count2 += 1   
            fileout.write(str(organism))
            fileout.write(str(framelabel))
            fileout.write('_')
            fileout.write(str(count2+1))
            fileout.write('_')
            fileout.write(str(poslist[count2]))
            fileout.write('\n')
            fileout.write(str(i))
            fileout.write('\n')
    # -r left at default runs orfsearch for all frames
    if args.readingframe == 0:
        orfsearch(capsequence, cutoff, 0)
        orfsearch(capsequence, cutoff, 1)
        orfsearch(capsequence, cutoff, 2)
        orfsearch(caprc, cutoff, 0)
        orfsearch(caprc, cutoff, 1)
        orfsearch(caprc, cutoff, 2)
    
    # allows users to select just 1 frame using -r option
    if args.readingframe == 1:
        orfsearch(capsequence, cutoff, 0)
    if args.readingframe == 2:
        orfsearch(capsequence, cutoff, 1)
    if args.readingframe == 3:
        orfsearch(capsequence, cutoff, 2)
    if args.readingframe == -1:
        orfsearch(caprc, cutoff, 0)
    if args.readingframe == -2:
        orfsearch(caprc, cutoff, 1)
    if args.readingframe == -3:
        orfsearch(caprc, cutoff, 2)
    
    #closing the output file - 1 file per input file
        fileout.close()

    #additional user statements printed to termninal to remind them of the cutoff and tell them the output file name
    print('with minimum ORF length of: ', args.cutoff)
    print('Printed to file: {0}'.format(filename))

exit()


# In[ ]:






