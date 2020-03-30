#!/usr/bin/env python3

import os
import json
import argparse

from convert_code_to_doc import *

def convertToMarkdownAndMerge(dir, config):

    for target, sources in config.items():

        targetExtension = os.path.splitext(target)[1]
        if not targetExtension == ".md":
            raise IOError("Markdown files expected as targets! Given target: {}".format(target))

        targetPath = os.path.join(dir, target)
        os.makedirs(os.path.dirname(targetPath), exist_ok=True)

        with open(targetPath, "w") as targetFile:
            for source in sources:
                fileExtension = os.path.splitext(source)[1]
                if fileExtension == ".md":
                    with open(os.path.join(dir, source), "r") as markdown:
                        targetFile.write(markdown.read())
                elif fileExtension == ".hh" or fileExtension == ".cc":
                    with open(os.path.join(dir, source), "r") as cppCode:
                        targetFile.write("\n\n" + transformCode(cppCode.read(), cppRules()) + "\n")
                else:
                    raise IOError("Unsupported or unknown file extension *{}".format(fileExtension))

def generateReadme(dir):
    config = None
    # look for .doc_config, if not found we pass
    try:
        configname = os.path.join(dir, ".doc_config")
        with open(configname, 'r') as configFile:
            config = json.load(configFile)
    except FileNotFoundError:
        pass
    if config is not None:
        convertToMarkdownAndMerge(dir, config)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="The folder to look for examples", default=".")
    args = vars(parser.parse_args())

    for path in os.listdir(args["directory"]):
        abspath = os.path.join(os.path.abspath(args["directory"]), path)
        if os.path.isdir(abspath):
            generateReadme(abspath)
