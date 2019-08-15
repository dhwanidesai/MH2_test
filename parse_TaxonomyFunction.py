import numpy as np
#import xarray as xr
import csv
from collections import defaultdict
import _pickle as cPickle
import bz2


import argparse, sys, textwrap

parser=argparse.ArgumentParser()

parser.add_argument('--taxafile', help = 'File mapping the reads to Taxa')
parser.add_argument('--taxafiletype', help = 'Source program for Taxa identification, e.g Kraken2, MEGAN, mmseqs2')
parser.add_argument('--funcfile', help = 'File mapping reads to functional categories')
parser.add_argument('--funcfiletype', help = 'Source program for Function identification, e.g mmseqs, diamond, MEGAN')
parser.add_argument('--m8file', help = 'BLAST or BLAST-like tab delimited output file mapping reads to Refseq or uniref IDs; required for Normalization to RPKG')

parser.add_argument('--multisample', help = textwrap.dedent('''For running multiple samples at a time, input a text file; overrides the above 5 arguments
    Format of the text file should contain 5 columns
    
    '''))

parser.add_argument('--MicrobeCensusReport', help = textwrap.dedent('''File containing Average genome size and genome equivalents info from MicrobeCensus program; 
Run the script runMicrobeCensus.py as follows to generate the report file:
python runMicrobeCensus.py --SampleDir test_reads/ --outfile microbeCensus-report.txt

Here --SampleDir is the directory containing the cleaned, stitched reads for the samples

If this argument is not supplied, the Normalization to RPKG will not be attempted
'''))
parser.add_argument('--unstratified', help='Output unstratified metabolic functions; must choose at most one of --stratified or --unstratified')
parser.add_argument('--stratified', help='Output metabolic functions stratified by taxa; must choose at most one of --stratified or --unstratified')


def main():
    args = parser.parse_args()
    
    multi = args.multisample
    strat = args.stratified
    unstrat = args.unstratified
    
    genlenf = bz2.BZ2File('RefGeneLength.pbz2', 'rb')
    genelendict = cPickle.load(genlenf)

    MCreport = args.MicrobeCensusReport
    GEdict = {}

    if (MCreport):
        GEdict = parseMicrobeCensusReport(MCreport)
        print (GEdict)
    else:
        print ("MicrobeCensus report not provided; RPKG will not be calculated")
    
    if not multi:
        taxadict = {}
        funcdict = {}

        taxafile = args.taxafile
        taxafiletype = args.taxafiletype
        funcfile = args.funcfile
        funcfiletype = args.funcfiletype
        m8file = args.m8file
        
        print ("Running single sample:", taxafile,taxafiletype,funcfile,funcfiletype)
        
        (taxadict,funcdict) = coreRun(taxafile,taxafiletype,funcfile,funcfiletype)
        
        print ("Total reads mapped to taxa: " + str(len(taxadict)))
        print ("Total reads mapped to functions: " + str(len(funcdict)))
        
        combinedDict = mergeTaxaFunc(taxadict,funcdict)
        print ("Reads mapped to either taxa OR functions: " + str(len(combinedDict)))
        filteredDict = dict(filter(lambda x: len(x[1]) == 2, combinedDict.items()))
                
        print ("Reads mapped to both taxa AND functions: " + str(len(filteredDict)))
 
        funchash = defaultdict(list)
        for key, value in sorted(filteredDict.items()):
            funchash[value[1]].append(key)
    
    
        funccount = mutate_dict(lambda x: len(x), funchash)
        
    else:
        print ("Running multiple samples from file:", multi)
        
        with open(multi, newline = '') as multif:                                                                                          
            multi_reader = csv.reader(multif, delimiter='\t')
            for line in multi_reader:
                print ("line:",line)
                    
                taxadict = {}
                funcdict = {}
                sampletag = line[0]
                taxafile = line[1]
                taxafiletype = line[2]
                funcfile = line[3]
                funcfiletype = line[4]
                m8file = line[5]

               # if (GEdict):
                #    GE = GEdict.get(sampletag)
                
                print ("Current sample:", taxafile,taxafiletype,funcfile,funcfiletype,m8file)
                (taxadict,funcdict,genedict ) = coreRun(taxafile,taxafiletype,funcfile,funcfiletype,m8file)
                print ("Total reads mapped to taxa: " + str(len(taxadict)))
                print ("Total reads mapped to functions: " + str(len(funcdict)))
        
                combinedDict = mergeTaxaFunc(taxadict,funcdict)
                print ("Reads mapped to either taxa OR functions: " + str(len(combinedDict)))
                filteredDict = dict(filter(lambda x: len(x[1]) == 2, combinedDict.items()))
                
                print ("Reads mapped to both taxa AND functions: " + str(len(filteredDict)))
 
                funchash = defaultdict(list)
                for key, value in sorted(filteredDict.items()):
                    funchash[value[1]].append(key)
            
            
                #first10pairs = {k: funchash[k] for k in list(funchash)[100:120]}
                #print("resultant dictionary : \n", first10pairs)
                print("Total unique functions: ", len(funchash))

                #FuncAvgGeneLengthdict = defaultdict(list)

                if (GEdict):
                    ge = GEdict.get(sampletag)
                else:
                    print ("MicrobeCensus reult not found; Will not normalize")
                
                for func,readarray in funchash.items():
                    FuncLengthArray = []
                    for read in readarray:
                        genes = genedict[read]
                        gene_lengths = list(map(lambda x: genelendict.get(x), genes))
                        list(np.float_(gene_lengths))
                        read_avg_gene_length = Average(gene_lengths)
                        #print (read,avg_gene_length)
                        #FuncAvgGeneLengthdict[read].append(avg_gene_length)
                        FuncLengthArray.append(read_avg_gene_length)

                    FuncAvgGeneLength = Average(FuncLengthArray)
                    #FuncAvgGeneLengthdict[func].append(FuncAvgGeneLength)

                    numreads = len(readarray)
                    if (ge):
                        lengthkb = FuncAvgGeneLength/1000
                        rpkg = numreads/lengthkb/ge

                    print ("Function: " + str(func) + " RPKG: " + str(rpkg))    
                #first10pairs = {k: FuncAvgGeneLengthdict[k] for k in list(FuncAvgGeneLengthdict)[100:120]}
                #print("resultant dictionary : \n", first10pairs)

                #print ()
                
                #if (GE):
                            
                #for func,readarray in funchash.items():
                   # for read in readarra
                
                            
                            
			
                #funccount = mutate_dict(lambda x: len(x), funchash)
                
                #if(unstrat):
                    
                    
                


#def unstratified():
    
    
#def stratified():
    

#def normalize_rpkg(dict_reads_mapped_to_func):
    

def parseMicrobeCensusReport(MCRfile):
    d = {}
    with open(MCRfile) as f:
        for line in f:
            fields = line.split("\t")
            S_tag = fields[0]

            if S_tag == "SampleTag":
                continue
            
            S_GE = float(fields[3])
            d[str(S_tag)] = S_GE    

    return d

def Average(lst):
    avg=0
#    print ("In average .... yay!")

    if(len(lst) > 0):
        avg = float(sum(np.float_(lst)) / len(lst))
            
    return avg

            
def mutate_dict(f,d):
    new_d = {}
    # apply a function to all elements of a dict
    # to generate a new dict
    for k, v in d.items():
        new_d[k] = f(v)
    
    return new_d

def coreRun(taxaf,taxaft,funcf,funcft,m8):
    
    taxadict = {}
    funcdict = {}

    print ("In Core ... >>>:",taxaft,funcft)
    
    if (taxaft == 'megan'):
        print ("Taxatype:Megan\n")
        taxadict=parseMeganTaxafile(taxaf)
    elif (taxaft == 'kraken2'):
        print ("Taxatype:Kraken2\n")
        taxadict=parseKraken2Taxafile(taxaf)
    else:
        print ("Taxa type not recognised\n")
        
        
    if (funcft == 'megan'):
        print ("Functype:Megan\n")
        funcdict = parseMeganFuncfile(funcf)
    elif (funcft == 'other'):
        print ("Functype:Other\n")
        funcdict = parseOtherFuncfile(funcf)
    else:
        print ("Func type not recognised\n")
    
    if (m8):
        print ("m8 file given; will Normalize to get RPKG")
        genedict = parseBlastm8(m8)
    
    return taxadict,funcdict,genedict
    

def mergeTaxaFunc(dict1,dict2):
    # Combine the values with same keys 
    result = defaultdict(list)
    # 
    for d in (dict1,dict2):
        for k,val in d.items():
            result[k].append(val)

            
    #first10pairs = {k: result[k] for k in list(result)[10:30]}
    #print("resultant dictionary : \n", first10pairs) 
    return result

def parseBlastm8(filename):
    d = defaultdict(list)
    print ("In m8 parse ... ")
    
    with open(filename) as f:
        for line in f:
            #(key, val) = line.split("\t")
            fields = line.split("\t")
            key = fields[0]
            val = fields[1] 
            #d[str(key)] = str(val)
            d[key].append(val)
    
    #first10pairs = {k: d[k] for k in list(d)[10:20]}
    #print("resultant dictionary : \n", first10pairs) 
    return d

    
    

def parseMeganFuncfile(filename):
    d = {}
    print ("In Megan func ... ")
    with open(filename) as f:
        for line in f:
            #fields = line.split("\t")
            (key, val) = line.split("\t")
            d[str(key)] = str(val)
    return d

def parseMeganTaxafile(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            (key, val) = line.split("\t")
            d[str(key)] = str(val)
    return d        

def parseKraken2Taxafile(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            fields = line.split("\t")
            classified = fields[0]
            key = fields[1]
            val = fields[2]    
            #print (fields[0])
            if (classified == 'C'):
                d[str(key)] = str(val)
    return d

if __name__ == "__main__":
    main();


    
