#!/usr/bin/env python3
import argparse
import sys
from converttomarkdown import convertToMarkdown

##########################################
# Convert source file into markdown file #
##########################################
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--source", help="Path to the source code file")
parser.add_argument("-t", "--target", help="Path to the target markdown file", default = "")
args = vars(parser.parse_args())

sourceFileName = args["source"]
targetFileName = args["target"]

if targetFileName == "":
    print( convertToMarkdown(sourceFileName) )

else:
    if targetFileName[-3:] != ".md":
        sys.stderr.write("Error: expected markdown file name as target\n")
        sys.exit(1)

    targetFile = open(targetFileName, 'w')
    targetFile.write( convertToMarkdown(sourceFileName) )
