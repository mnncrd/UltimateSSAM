"""Ultimate Secondary Structure Assignment Method.

This program assigns secondary structures to a sequence.

    Typical usage example:

    $ python ssam.py filename -o output_file
"""

import argparse

def main():

    """Main function"""

    #Get the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="The file to assign a secondary structure to", type=str)
    parser.add_argument("-o", type=str)
    parser.add_argument("-hy", "--hydrogen", action="store_true", help="add hydrogen atoms")
    args = parser.parse_args()
