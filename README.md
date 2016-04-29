python NCBITK.py -h
usage: NCBITK.py [-h] [-W] [-i INPUT_FILE] [-l FROM_LIST] local_mirror

Sync with NCBI's database, give the files useful names,and organize them in a
sane way.

positional arguments:
  local_mirror          Your local directory to save fastas to, e.g.
                        "bacteria"

optional arguments:
  -h, --help            show this help message and exit
  -W, --no_wget         Don't fetch assembly_summary.txt
  -i INPUT_FILE, --input_file INPUT_FILE
                        Input file containing directories to sync with.
  -l FROM_LIST, --from_list FROM_LIST
                        Comma separated list of directories to be downloaded
