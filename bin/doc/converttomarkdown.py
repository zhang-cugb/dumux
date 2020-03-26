#!/usr/bin/env python3
import sys

#####################################################
# Function to create markdown text from source code #
#####################################################
def convertToMarkdown(sourceFileName):

    # check admissible file extensions
    if sourceFileName[-3:] == ".hh":
        isHeader = True
    elif sourceFileName[-3:] == ".cc":
        isHeader = False
    else:
        sys.stderr.write("Error: supported file formats: \".hh\" and \".cc\"\n")
        sys.exit(1)

    result = ""
    startParsing = False
    insideCodeBlock = False
    insideLargeCodeBlock = False

    for line in open(sourceFileName, 'r').read().splitlines():

        # keep reading until the actual start of the code
        if startParsing == False:
            if isHeader == True and ("#define" in line):
                startParsing = True
            elif isHeader == False and ("***/" in line):
                startParsing = True

        # parse lines of source code
        else:
            # begin a large code block
            if insideLargeCodeBlock == False and ("<codeblock>" in line):
                if insideCodeBlock == True:
                    result += "```\n"
                    insideCodeBlock = False
                insideLargeCodeBlock = True
                result += "\n```cpp\n"

            # inside a large code block
            elif insideLargeCodeBlock == True:
                if ("</codeblock>" in line):
                    insideLargeCodeBlock = False
                    result += "\n```\n"
                else:
                    result += line + "\n"

            # Here, all comments are passed as markdown & the rest as code
            else:
                # possibly descriptive text for the markdown file
                if "//" in line:
                    # make sure this is not an in-line comment inside code
                    line = line.split("//")

                    # to be added as markdown
                    if len(line) == 1 \
                       or (len(line) > 1 and all(c in " \t" for c in line[0])) \
                       or (len(line) > 1 and line[0] == ""):

                        # end previous code block
                        if insideCodeBlock == True:
                            result += "```\n"
                            result += "\n" + line[ len(line)-1 ] + "\n"
                            insideCodeBlock = False
                        else:
                            result += line[ len(line)-1 ] + "\n"

                    # in-line comment within code
                    elif len(line) == 2:
                        if insideCodeBlock == False:
                            result += "\n```cpp\n"
                            insideCodeBlock = True
                        result += line[0] + "//" + line[1] + "\n"

                # this is part of code
                else:
                    if insideCodeBlock == False and line != "":
                        result += "\n```cpp\n"
                        insideCodeBlock = True
                    result += line + "\n"
    return result
