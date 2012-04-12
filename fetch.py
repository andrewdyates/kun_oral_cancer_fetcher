#!/usr/bin/python
"""
Write a list of cleaned gene expression data to compressed files.

Input list as plain text in STDIN, one line per GSE ID.

EXAMPLE USE:

$ cat my-gse-list.txt | python fetch.py 2> mylog.txt
"""
import geo_api
import sys
import os
import gzip

def main():
  # Fetch list of GSEs from STDIN
  gse_list = [line.strip() for line in filter(None, sys.stdin.readlines())]
  fetch_gse_list(gse_list)


def fetch_gse_list(gse_list, out_dir=""):
  """Fetch a list of expression GSEs, write to compressed files.

  Args:
    gse_list: [str] of GSE\d+ IDs of expression-only GEO studies.

  Assume that no GPLs or substudies must be specified.
  """
  for gse_id in gse_list:
    g = geo_api.GSE(gse_id)
    # if this is a super study, select the child that is an eQTL study
    if g.type == "SUPER":
      for gg in g.substudies.values():
        if gg.type == "eQTL":
          break
      g = gg
    if g.type != "eQTL":
      print "Could not find eQTL study for %s. Skipping..." % gse_id
      continue
    filt = geo_api.EQTLFilter(g)
    filename = "%s.clean.tab.gzip" % gse_id
    filepath = os.path.join(out_dir, filename)
    fp = gzip.open(filepath, "wb")
    print "Writing %s..." % filepath
    nn = None
    for row in filt.get_rows():
      # First 5 columns: 'ID_REF', 'GENE_SYMBOL', 'NUM_VALUES', 'MEAN', 'STD'
      #   remove all but 'GENE_SYMBOL' column which should be made first column
      row = [row[1]] + row[5:]

      # Verify column alignment
      n = len(row)
      if nn is not None and n != nn:
        raise Exception, "Columns not aligned."
      nn = n

      # Replace None with empty string.
      row = map(lambda x: x if x!="None" else "", row)

      # Print row. 
      fp.write("\t".join(row)); fp.write("\n")
    fp.close()
    
    
if __name__ == "__main__":
  main()
