#!/usr/bin/python
import pandas as pd

def main():
    pathway_file = pd.read_csv("c2.cp.v6.2.entrez.gmt.csv")
    print(pathway_file)

if __name__ == '__main__':
    main()