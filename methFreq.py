#! /usr/bin/python
# Python 2.7.5, requires Biopython.

'''
Created on 07/06/2014
@author: Adam_Taranto
'''
import argparse;
import getopt;
import sys; 
import re; 
from Bio import SeqIO; 
from Bio.Seq import Seq;
#from Bio.Seq import MutableSeq;
from Bio.SeqRecord import SeqRecord;
#from Bio.Alphabet import IUPAC;

def main(filename, headerRow, decimalPlaces):

    # Set variables
    outFunc = writeTab
    headFunc = writeTab
    formatString = "{0:." + str(decimalPlaces) + "f}"
   
    # Do the work
    with open(filename, "rU") as handle:

      if headerRow:
          headFunc(("SeqID",
                    "Seq_len",
                    "GC_CONTENT",
                    "N_Count",
                    "CpG_OBS",
                    "CpG_EXP",
                    "CpG_OE"
                    "CpHpG_OBS",
                    "CpHpG_EXP",
                    "CpHpG_OE",
                    "CpHpH_OBS",
                    "CpHpH_EXP",
                    "CpHpH_OE"
                    ))
      
      # For each sequence
      for record in SeqIO.parse(handle, "fasta"):

          # Remove lowercase characters
          record = record.upper()

          # Convert mRNA to DNA
          if "U" in record.seq :
              toDNA = record.seq.back_transcribe()
              record = SeqRecord(toDNA, id=record.id, name=record.name)

          # Initialise counters!
          baseTotals = {'A':0.0, 'T':0.0, 'C':0.0, 'G':0.0,}

          pairTotals = {'AA':0.0,'AT':0.0,'AC':0.0,'AG':0.0, 
                        'TA':0.0,'TT':0.0,'TC':0.0,'TG':0.0, 
                        'CA':0.0,'CT':0.0,'CC':0.0,'CG':0.0, 
                        'GA':0.0,'GT':0.0,'GC':0.0,'GG':0.0,}
          
          tripletTotals = {'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0,
                           'ACA':0.0,'ACC':0.0,'ACG':0.0,'ACT':0.0,
                           'AGA':0.0,'AGC':0.0,'AGG':0.0,'AGT':0.0,
                           'ATA':0.0,'ATC':0.0,'ATG':0.0,'ATT':0.0,
                           'CAA':0.0,'CAC':0.0,'CAG':0.0,'CAT':0.0,
                           'CCA':0.0,'CCC':0.0,'CCG':0.0,'CCT':0.0,
                           'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,
                           'CTA':0.0,'CTC':0.0,'CTG':0.0,'CTT':0.0,
                           'GAA':0.0,'GAC':0.0,'GAG':0.0,'GAT':0.0,
                           'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,
                           'GGA':0.0,'GGC':0.0,'GGG':0.0,'GGT':0.0,
                           'GTA':0.0,'GTC':0.0,'GTG':0.0,'GTT':0.0,
                           'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,
                           'TCA':0.0,'TCC':0.0,'TCG':0.0,'TCT':0.0,
                           'TGA':0.0,'TGC':0.0,'TGG':0.0,'TGT':0.0,
                           'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0,}

          Sec2lastbase = 'N'
          lastbase = 'N'
          baseidx = 0
          seqlen = len(record.seq)

          # For each base in sequence
          for base in record.seq:
              # Sum the new triplet
              if base != 'N':
                # If current base is not N, then count current base
                baseTotals[base] += 1
                if lastbase != 'N':
                  #If lastbase was also not an N then count the current pair
                  pairTotals[lastbase+base] += 1
                  if Sec2lastbase != 'N':
                    #If no Ns in triplet count current triplet
                    tripletTotals[Sec2lastbase+lastbase+base] += 1

              # End of gene?
              if baseidx == (seqlen - 1):

                  # Calculate stats
                  gcContent = None
                  CpG_OBS = None
                  CpG_EXP = None
                  CpG_OE = None
                  CpHpG_OBS = None
                  CpHpG_EXP = None
                  CpHpG_OE = None
                  CpHpH_OBS = None
                  CpHpH_EXP = None
                  CpHpH_OE = None
                  
                  # Ncount
                  # Check calc with pair totals
                  nCount = len(record.seq) - sum(baseTotals.values())
                  ##nCount = len(record.seq) - sum(pairTotals.values())

                  if seqlen > 0:
                      probA = baseTotals['A'] / seqlen
                      probT = baseTotals['T'] / seqlen
                      probG = baseTotals['G'] / seqlen
                      probC = baseTotals['C'] / seqlen
                      probH = (1 - baseTotals['G']) / seqlen

                  # Only need to validate that you are not going to divide by zero.

                  # CpG_OBS
                  if seqlen > 0:
                      CpG_OBS_Num = pairTotals['CG']
                      CpG_OBS = formatString.format(CpG_OBS_Num)

                  # CpG_EXP
                  if seqlen > 0:
                      CpG_EXP_Num = (baseTotals['C'] * baseTotals['G']) / seqlen
                      CpG_EXP = formatString.format(CpG_EXP_Num)

                      # CpG_OE
                      if CpG_EXP_Num > 0:
                          CpG_OE_Num = CpG_OBS_Num / CpG_EXP_Num
                          CpG_OE = formatString.format(CpG_OE_Num)

                  ## CpHpG
                  if seqlen > 0:
                      CpHpG_OBS_Num = tripletTotals['CAG'] + tripletTotals['CCG'] + tripletTotals['CTG']
                      ##CpHpG_EXP_Num_NEG = tripletTotals['CTG'] + tripletTotals['CGG'] + tripletTotals['CAG']
                      CpHpG_OBS = formatString.format(CpHpG_OBS_Num)

                  # CpHpG_EXP
                  if seqlen > 0:
                      CpHpG_EXP_Num = (probC * probH * probG) * seqlen
                      CpHpG_EXP = formatString.format(CpHpG_EXP_Num)

                      # CpHpG_OE
                      if CpHpG_EXP_Num > 0:
                          CpHpG_OE_Num = CpHpG_OBS_Num / CpHpG_EXP_Num
                          CpHpG_OE = formatString.format(CpHpG_OE_Num)

                  ## CpHpH
                  if seqlen > 0:
                      
                      CpHpH_OBS_Num = (tripletTotals['CAA'] + tripletTotals['CCA'] + tripletTotals['CTA'] + 
                                       tripletTotals['CAT'] + tripletTotals['CCT'] + tripletTotals['CTT'] + 
                                       tripletTotals['CAC'] + tripletTotals['CCC'] + tripletTotals['CTC'])  
                      ##CpHpH_OBS_Num_NEG = (tripletTotals['GAG'] + tripletTotals['GGG'] + tripletTotals['GTG'] + 
                      ##                     tripletTotals['AAG'] + tripletTotals['AGG'] + tripletTotals['ATG'] + 
                      ##                     tripletTotals['TAG'] + tripletTotals['TGG'] + tripletTotals['TTG'])
                      CpHpH_OBS = formatString.format(CpHpH_OBS_Num)

                  # CpHpH_EXP
                  if seqlen > 0:
                      CpHpH_EXP_Num = (probC * probH * probH) * seqlen
                      CpHpH_EXP = formatString.format(CpHpH_EXP_Num)

                      # CpHpH_OE
                      if CpHpH_EXP_Num > 0:
                          CpHpH_OE_Num = CpHpH_OBS_Num / CpHpH_EXP_Num
                          CpHpH_OE = formatString.format(CpHpH_OE_Num)

                  # GC content    
                  ACGT = baseTotals['A'] + baseTotals['T'] + baseTotals['C'] + baseTotals['G']
                  if ACGT > 0:
                      gcContent = formatString.format((baseTotals['C'] + baseTotals['G']) / ACGT)
                  
                  # Print the results
                  outFunc((record.id,                 # seq id
                            str(seqlen),              # sequence length
                            str(gcContent),           # GC content
                            str(nCount),              # total N's in window
                            str(CpG_OBS),
                            str(CpG_EXP),
                            str(CpG_OE),
                            str(CpHpG_OBS),
                            str(CpHpG_EXP),
                            str(CpHpG_OE),
                            str(CpHpH_OBS),
                            str(CpHpH_EXP),
                            str(CpHpH_OE),
                            ))

              # Update counters and trackers
              baseidx += 1
              Sec2lastbase = lastbase
              lastbase = base

          # Loop to next base
      # next sequence
    handle.close()

## Output format writer functions ##
def writeTab(record):
    '''Writes a record in Tab-delimited format'''
    delimiter = "\t"
    print(delimiter.join(record))

if __name__=='__main__':
  ### Argument handling
  arg_parser = argparse.ArgumentParser(description='Script description');
  arg_parser.add_argument("filename", help="A fasta file containing DNA coding sequences");
  arg_parser.add_argument("-H", "--header", type=bool, default=True, help="Print header row on tab delimited output");
  arg_parser.add_argument("-d", "--decimal", type=int, default=3, help="Format values to x decimal places");
  args = arg_parser.parse_args();

  ### Variable definitions/declarations
  filename = args.filename;
  headerRow = args.header;
  decimalPlaces = args.decimal;
 
  ## Pass variables to main script
  main(filename, headerRow, decimalPlaces);
