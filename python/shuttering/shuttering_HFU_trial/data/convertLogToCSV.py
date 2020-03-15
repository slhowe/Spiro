#!/usr/bin/env python
import csv

# Got a log file to convert to csv format
#
# Input file format is European numbers separated by spaces
# eg 0,09 0,8 1,2
# Output file format is  English numbers separated by commas
# eg 0.09,0.8,1.2

def main():
    inputFileName = "DIFFUSTIK.log"
    outputFileName = "DIFFUSTIK.csv"

    with open(inputFileName, 'r') as f:
        with open(outputFileName, 'w') as outputFileName:

            header = f.readline()

            # This is the first line in the output file
            header = ["Flow","PM","CO","CH4"]

            writer = csv.writer(outputFileName, delimiter=",")
            writer.writerow(header)

            for line in f:
                line = line.replace(',', '.')
                line  = line.split()
                reformatted_line = []
                for num in line:
                    num = '  ' + num
                    reformatted_line.append(num)
                writer.writerow(reformatted_line)
        outputFileName.close
    f.close

if __name__ == "__main__":
    main()
