import os
import sys
import argparse
from pdfrw import PdfReader, PdfWriter 




def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="File to inspect for profiler information. Currently, only PDF is supported.")
    args = parser.parse_args()

    if not os.path.isfile(args.file):
        sys.stderr.write("Path seems to not exist or is not a regular file!")


    trailer = PdfReader(args.file)
    for key in trailer.Info.keys():
        sys.stdout.write('Key "{}"\n'.format(key))
        sys.stdout.write("{}\n\n".format(trailer.Info[key]))

#end main()

if __name__ == "__main__":
    main()