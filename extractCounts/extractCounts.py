#!/usr/bin/python
import argparse, straw, sys
import pandas as pd

if __name__ == "__main__":

   ## Initialize argument parser
   parser = argparse.ArgumentParser(
      description = "extractCounts.py uses straw to extract counts from a list of Hi-C files at loop positions in bedpe format."
   )

   ## Add arguments
   parser.add_argument('-l', '--loopFile', required=True, help='BEDPE file of loop locations')
   parser.add_argument('-f', '--hicFiles', required=True, help='List of Hi-C files', nargs='+')
   parser.add_argument('-o', '--outFile', required=False, help='File to write output; default is loopsOut.txt', default="loopsOut.txt")
   parser.add_argument('-n', '--norm', required=False, help='Normalization <NONE/VC/VC_SQRT/KR>; default is NONE', default='NONE')

   ## Parse command-line arguments
   args = parser.parse_args()

   ## List Hi-C files from which counts should be extracted
   hicFiles = args.hicFiles

   ## Read in loops files & convert all columns to strings
   loops = pd.read_table(args.loopFile).astype(str)

   ## Concatenate pos1 and pos2 vectors (also remove "chr" from chromosome name)
   pos1 = loops.iloc[:,0].replace('chr', '', regex=True) + ":" +  loops.iloc[:,1] + ":" +  loops.iloc[:,1]
   pos2 = loops.iloc[:,3].replace('chr', '', regex=True) + ":" +  loops.iloc[:,4] + ":" +  loops.iloc[:,4]

   ## Resolution for each loop
   res = loops.iloc[:,2].astype(int) - loops.iloc[:,1].astype(int)

   ## Add counts to loops file for each Hi-C file
   for h in hicFiles:
      for i in range(len(loops)):
         try:
               loops.loc[i, h] = straw.straw(args.norm, h, pos1[i], pos2[i], 'BP', res[i])[2][0]
         except:
               loops.loc[i, h] = 0
               print('row {} in file {} was imputed with 0'.format(i,h), file=sys.stderr)

   # Write as tab-separate file
   loops.to_csv(args.outFile, sep='\t')