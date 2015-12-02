
from __future__ import print_function

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import argparse
import gzip
import re
import os
import math

#plt.interactive(True)

def main():

    parser=argparse.ArgumentParser(description='signal pile up plot for bedGraph/wiggle tracks',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='signalTrack', type=str, required=True, help='input bed/wiggle track to be used for the pile up')
    parser.add_argument('-e', '--element', dest='elementTrack', type=str, required=True, help='element bed track for pile up')
    parser.add_argument('-g', '--genome', dest='genome', type=str, required=True, help='genome chromosome length file')
    parser.add_argument('-w', '--width', dest='pileUpWidth', type=int, required=False, default=500, help='width of pile up plot (in BP)')
    parser.add_argument('-b', '--num_bins', dest='num_bins', type=int, required=False, default=500, help='number of num_bins for pile up plot width')
    parser.add_argument('-ya','--yaxisrange', dest='yaxisrange', type=int, nargs='+', required=False, default=[0,0], help='y-axis for aggregrate plot')
    parser.add_argument('-sc', '--signalcolumn', dest='signalColumn', type=int, required=False, default=3, help='signal column number in bed/bedgraph/wig/tsv file')
    parser.add_argument('-sm', '--sortMatrix',dest='sortMatrix', action='store_true', help='sort pile up matrix by total row signal, descending')
    parser.add_argument('-mp', '--midpointMode',dest='midpointMode', action='store_true', help='use midpoint of element, otherwise will assume TSS')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    
    args=parser.parse_args()
    
    signalTrack=args.signalTrack
    elementTrack=args.elementTrack
    genome=args.genome
    pileUpWidth=args.pileUpWidth
    num_bins=args.num_bins
    yaxisrange=args.yaxisrange
    signalColumn=args.signalColumn
    sortMatrix=args.sortMatrix
    midpointMode=args.midpointMode
    verbose=args.verbose
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    verboseprint("")
    verboseprint("signalTrack",signalTrack,sep="\t")
    verboseprint("elementTrack",elementTrack,sep="\t")
    verboseprint("genome",genome,sep="\t")
    verboseprint("pileUpWidth",pileUpWidth,sep="\t")
    verboseprint("num_bins",num_bins,sep="\t")
    verboseprint("yaxisrange",yaxisrange,sep="\t")
    verboseprint("signalColumn",signalColumn,sep="\t")
    verboseprint("sortMatrix",sortMatrix,sep="\t")
    verboseprint("midpointMode",midpointMode,sep="\t")
    verboseprint("verbose",verbose,sep="\t")
    verboseprint("")
    
    verboseprint("validating element vector")
    element_size=ensure_sorted(elementTrack)
    verboseprint("\t",element_size," elements\n",sep="")
    
    verboseprint("validating signal vector")
    signal_size=ensure_sorted(signalTrack)
    verboseprint("\t",signal_size," signals\n",sep="")
    
    # standardize names
    signalTrack_name=os.path.basename(signalTrack)
    signalTrack_name=re.sub(".gz", "", signalTrack_name)
    elementTrack_name=os.path.basename(elementTrack)
    elementTrack_name=re.sub(".gz", "", elementTrack_name)
    pileUpName=signalTrack_name+"__"+elementTrack_name
    
    # read through the genome file
    chr_dict,genomeSize=process_genome(genome)
    
    # process element file - get element list
    verboseprint("processing elements")
    elements=process_elements(elementTrack,midpointMode,chr_dict,pileUpWidth,verboseprint)
    element_size=len(elements)
    verboseprint("\tdone\n")
    
    # read through the signal tracks
    verboseprint("processing signal vector")
    if signalTrack.endswith('.gz'):
        f1=gzip.open(signalTrack,"r") 
    else:
        f1=open(signalTrack,"r") 
    
    bin_size=int(max(math.ceil(((pileUpWidth*2)+1)/num_bins),1))
    if((bin_size % 2) != 0): # force bin_size to be even
        bin_size = bin_size + 1 
    num_bins=int(math.ceil((pileUpWidth*2)/bin_size))
    if((num_bins % 2) == 0): # force num_bins to be odd
        num_bins=num_bins + 1 
    pileUpWidth=((num_bins*bin_size)/2)
    
    elementSignal=np.zeros(num_bins,dtype=np.float32)
    elementCount=np.zeros(num_bins,dtype=np.int)
    counts=np.zeros(num_bins,dtype=np.int)
    pileUpMatrix=np.zeros((element_size,num_bins),dtype=np.float32)
    pileUpMatrix.fill(np.nan) # fill with NaN
    
    verboseprint("\tpileUpWidth\t",-pileUpWidth,"bp - +",pileUpWidth,"bp",sep="")
    verboseprint("\tnBins\t",num_bins,sep="")
    verboseprint("\tbin_size\t",bin_size,"bp",sep="")
    
    lastElementIndex=0
    signalIndex=0
    signalBucketSize=max(1,int(signal_size/1000))
    elementBucketSize=max(1,int(element_size/1000))
    
    elementIndex=0
    tmpElement=elements[elementIndex]
    where=None
    init=0
    
    # iterate through signal file
    while True:
        line = f1.readline()
        
        # if EOF, break
        if not line: 
            break 
      
        line = line.rstrip("\n")
        
        if line.startswith("track") or line.startswith("#"): # skip track def, or comment lines
            continue
        
        tmpArr=line.split("\t")
        
        # assuming bedGraph/wig format (chr,start,end,score)
        chromosome=tmpArr[0]
        start=int(tmpArr[1])
        end=int(tmpArr[2])-1 # assume half open interval (a,b]
        signal=float(tmpArr[signalColumn])
        
        #print(where,tmpArr,tmpElement,lastElementIndex,elementIndex)
            
        if((lastElementIndex != -1) and (lastElementIndex != elementIndex)):
                    
            # reset where to f1 
            print("RESET WHERE",where,f1.tell())
            where = f1.tell()

            print("NEW ELEMENT",where,lastElementIndex,elementIndex,elements[lastElementIndex],elements[elementIndex])
            
            nan_idx = elementSignal==0
            elementSignal[nan_idx]=np.nan
            elementCount_mask = elementCount > 0
            counts = counts + elementCount_mask
            
            with np.errstate(invalid='ignore'):
                meanSignal=np.divide(elementSignal,elementCount)
                
            pileUpMatrix[lastElementIndex,:]=meanSignal
            print("logging pileupmatrix",lastElementIndex,"=",elementSignal)
            
            # reset buffer
            elementSignal=np.zeros(num_bins,dtype=np.float32)
            elementCount=np.zeros(num_bins,dtype=np.int)
        
        
        overlap=0
        if(chromosome == tmpElement[0]):
            overlap=is_overlap([tmpElement[1],tmpElement[3]],[start,end])
            
            if(overlap > 0):   # if we have an overlap
                
                print("\toverlap",init,overlap,lastElementIndex,elementIndex)
                
                if(init == 0):
                    # reset where to f1 
                    print("RESET WHERE",where,f1.tell())
                    where = f1.tell()
                init=1
                
                # get offset into element, handle bounds
                startOffset=(start-tmpElement[1])
                endOffset=(end-tmpElement[1])
                startIdx=max(0,(startOffset/bin_size))
                endIdx=min((num_bins-1),(endOffset/bin_size))
                print("\tOFFSETS",startOffset,startIdx,endOffset,endIdx)
                
                # aggregrate signal data
                x = np.arange(startIdx,endIdx+1)
                print("\tOVERLAP",num_bins,bin_size,(pileUpWidth*2),start,end,tmpElement,startOffset,endOffset,x)
                elementSignal[x] = elementSignal[x] + signal
                print("\t\tlogging",signal,np.nansum(elementSignal[x]))
                elementCount[x] = elementCount[x] + 1
                print(elementCount)
                    
        
        if(((signalIndex % signalBucketSize) == 0) or (signalIndex == (signal_size-1)) or ((lastElementIndex != elementIndex) and ((elementIndex % elementBucketSize) == 0))):
            element_pc=((float(elementIndex+1)/float((element_size)))*100)
            signal_pc=((float(signalIndex+1)/float((signal_size)))*100)
            verboseprint("\r\telement ["+str(tmpElement[0])+"] "+str(elementIndex)+"/"+str(element_size-1)+" ["+str("{0:.2f}".format(element_pc))+"%] | signal ["+str(tmpArr[0])+"] "+str(signalIndex)+"/"+str(signal_size-1)+" ["+str("{0:.2f}".format(signal_pc))+"%] complete ... ",end="\r")
            if verbose: sys.stdout.flush()
        
        lastElementIndex=elementIndex
        
        while( ( (chromosome > tmpElement[0]) or ((chromosome == tmpElement[0]) and  (start > tmpElement[3])) ) and (elementIndex < (element_size-1))):
            
            init=0
            
            elementIndex = elementIndex+1
            print("next element",lastElementIndex,"->",elementIndex,tmpElement,chromosome,start,end,signalIndex)
            tmpElement=elements[elementIndex]
            
            element_overlap=0
            if((lastElementIndex != -1) and (elements[lastElementIndex][0] == elements[elementIndex][0])):
                element_overlap=is_overlap([elements[lastElementIndex][1],elements[lastElementIndex][3]],[elements[elementIndex][1],elements[elementIndex][3]])
            
            print(elements[lastElementIndex],elements[elementIndex],element_overlap,where)
            
            if(element_overlap != 0):
                if where is not None:
                    print("\nabout to seek",lastElementIndex,elementIndex,element_overlap,where,signalIndex)
                    f1.seek(where)
                    print("done seek\n")
                    break
        
        signalIndex = signalIndex + 1
        
    f1.close()
    
    # clean up buffer for last element
    if(lastElementIndex != -1):
        
        nan_idx = elementSignal==0
        elementSignal[nan_idx]=np.nan
        elementCount_mask = elementCount > 0
        counts = counts + elementCount_mask
        
        with np.errstate(invalid='ignore'):
            meanSignal=np.divide(elementSignal,elementCount)
        
        #print("logging pileupmatrix",lastElementIndex,"=",elementSignal)
        pileUpMatrix[lastElementIndex,:]=meanSignal
            
        element_overlap=0
        if((lastElementIndex != -1) and (elements[lastElementIndex][0] == elements[elementIndex][0])):
            element_overlap=is_overlap([elements[lastElementIndex][1],elements[lastElementIndex][3]],[elements[elementIndex][1],elements[elementIndex][3]])
        
    # calculate aggregrate signal as sum
    aggregrate=np.nansum(pileUpMatrix,axis=0)
    aggregrate=aggregrate/counts
    #print(aggregrate,counts)
    
    np.set_printoptions(threshold='nan')
    #print(aggregrate)
    
    if((len(yaxisrange) != 2) or (yaxisrange[0] == yaxisrange[1])):
        yaxisrange=[0,max(aggregrate)]

    # plot the aggregrate signal
    bin_starts = np.arange(-pileUpWidth,pileUpWidth,bin_size)
    bin_midpoints = np.arange(-pileUpWidth+(bin_size/2),pileUpWidth+(bin_size/2),bin_size)
    #print(len(bin_starts),len(bin_midpoints))
    #print(bin_starts,bin_midpoints)
    
    # open output file
    agg=gzip.open(pileUpName+'.aggregrate.vector.gz',"wb")
    for i,bin in enumerate(bin_midpoints):
        print(bin_starts[i],bin_starts[i]+bin_size,bin,aggregrate[i],sep="\t",file=agg)
    print("\n",sep="",end="",file=agg)
    agg.close()
    
    plt.plot(bin_midpoints,aggregrate)
    plt.xlabel('genome distance from element-anchor')
    plt.ylabel('aggregrate signal')
    plt.title(signalTrack_name+"\n"+elementTrack_name)
    plt.grid(True)
    axisbounds=[-pileUpWidth,pileUpWidth]+yaxisrange
    plt.axis(axisbounds)
    plt.savefig(pileUpName+'.aggregrate.vector.png')

    # row sums of matrix
    rowsums=np.nansum(pileUpMatrix,axis=1)
    
    # sort by rowsums, get top 25% or 5000 idx, descending order
    idx = np.arange(0,len(rowsums))
    if(sortMatrix == 1):
        idx = rowsums.argsort()[::-1][:max(5000,int(pileUpMatrix.shape[0]*.25))]
    
    # open output file
    out_fh=gzip.open(pileUpName+'.matrix.gz',"wb")
    
    print("# signalTrack",signalTrack,sep="\t",file=out_fh)
    print("# elementTrack",elementTrack,sep="\t",file=out_fh)
    print("# genome",genome,sep="\t",file=out_fh)
    print("# pileUpWidth",pileUpWidth,sep="\t",file=out_fh)
    print("# num_bins",num_bins,sep="\t",file=out_fh)
    print("# yaxisrange",yaxisrange,sep="\t",file=out_fh)
    print("# signalColumn",signalColumn,sep="\t",file=out_fh)
    print("# sortMatrix",sortMatrix,sep="\t",file=out_fh)
    print("# midpointMode",midpointMode,sep="\t",file=out_fh)
    print("# verbose",verbose,sep="\t",file=out_fh)
    
    for x in xrange(-pileUpWidth,pileUpWidth,bin_size):
        print("\t","pos|",x,"__",x+(bin_size/2),"__",x+bin_size,sep="",end="",file=out_fh)
    print("\n",sep="",end="",file=out_fh)
    
    for i,v in enumerate(pileUpMatrix[idx,:]):
        if((sortMatrix == 1) and (rowsums[idx[i]] == 0)): # remove blank rows, if sort mode enabled
            continue
        print("element_",elements[idx[i]][5],"\t","\t".join(map(str,v)),sep="",file=out_fh) 
    
    out_fh.close()
                    
    verboseprint("")
    verboseprint("")
    
def ensure_sorted(file):
    """ensure file is sorted by chr,start,end - return nLines
    """
    
    if file.endswith('.gz'):
        fh=gzip.open(file,'r')
    else:
        fh=open(file,"r")
    
    lastBED3=[]
    nData=0;
    for i,line in enumerate(fh):
        line=line.rstrip("\n")
        
        if line.startswith("track") or line.startswith("#"):
            continue
        
        tmpArr=line.split("\t")
        chr=tmpArr[0]
        start=int(tmpArr[1])
        end=int(tmpArr[2])
        
        if len(lastBED3) == 0:
            lastBED3=[chr,start,end]
        
        if chr < lastBED3[0]:
            sys.exit('\nsort disorder  (chr)! line# '+str(i)+' | '+chr+' < '+str(lastBED3[0])+'!\n\t'+file+'\n')
        if (chr == lastBED3[0]) and (start < lastBED3[1]):
            sys.exit('\nsort disorder (start)! line# '+str(i)+' | '+str(start)+' < '+str(lastBED3[1])+'!\n\t'+file+'\n')
        if (chr == lastBED3[0]) and (start == lastBED3[1]) and (end < lastBED3[2]):
            sys.exit('\nsort disorder (end)! line# '+str(i)+' | '+str(end)+' < '+str(lastBED3[2])+'!\n\t'+file+'\n')
        
        lastBED3=[chr,start,end]
        nData=nData+1
        
    return(nData)

def process_genome(genome):
    """process the genome file, build chr length dict, and genome size
    """

    genomeSize=0
    chr_dict={}
    
    g_fh=open(genome,"r")
    for line in g_fh:
        line=line.rstrip("\n")
        
        if line.startswith("track") or line.startswith("#"):
            continue
        
        tmpArr=line.split("\t")
        
        start=1
        if(len(tmpArr) == 3):
            start=int(tmpArr[1])
            end=int(tmpArr[2])
        
        if(len(tmpArr) == 2):
            end=int(tmpArr[1])
        
        # store chromosome start,end in dict
        chr_dict[tmpArr[0]]=[start,end]
        genomeSize = genomeSize + end
    g_fh.close()
    
    return(chr_dict,genomeSize)

def process_elements(elementTrack,midpointMode,chr_dict,pileUpWidth,verboseprint):
    """process the element file, return element list
    """

    if elementTrack.endswith('.gz'):
        et_fh=gzip.open(elementTrack,'r')
    else:
        et_fh=open(elementTrack,"r")
        
    skipped_elements=0
    elements=[]
    for i,line in enumerate(et_fh):
        line=line.rstrip("\n")
        
        if line.startswith("track") or line.startswith("#"):
            continue

        tmpArr=line.split("\t")
        chromosome=tmpArr[0]
        strand="+"
        if(len(tmpArr) >= 6):
            strand=tmpArr[5]
        start=int(tmpArr[1])
        end=int(tmpArr[2])
        name=tmpArr[3]
        midpoint=int((start+end)/2) # always round up
        
        anchor=midpoint
        if not midpointMode: # if non-midpoint mode, use TSS 
            if(strand == "+"):
                anchor = start
            else:
                anchor = end
            
        # get chr bounds
        if chromosome not in chr_dict:
            skipped_elements=skipped_elements+1
            continue
            
        # skip chrs not in 
        chrStart,chrEnd=chr_dict[chromosome]
        
        # stretch out anchor by pileUpWidth, limit by chr bound
        newStart=max(chrStart,(anchor-pileUpWidth))
        newEnd=min(chrEnd,(anchor+pileUpWidth))
        
        # append element obj to list
        elements.append([chromosome,newStart,anchor,newEnd,strand,name])
        
    et_fh.close()
    
    verboseprint("\tskipped ",skipped_elements," elements",sep="")
    
    return(elements)

def is_overlap(a, b):
    """test to for overlap between two intervals.
    """
    
    if(a[0] > a[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(a[0])+' > end '+str(a[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    if(b[0] > b[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(b[0])+' > end '+str(b[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
        
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

if __name__=="__main__":
    main()

   