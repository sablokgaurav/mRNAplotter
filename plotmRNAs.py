def plotmRNAs(inputgff, outfile):
    """
    Author Gaurav Sablok
    Universitat Potsdam
    Date: 2024-4-16
    plotting the mRNA for the genome annotation and making
    the graph layout of the mRNAs using the graphql API. 
    @inputfile = file aligned to the genome using the protein hints
    @outputfile = file to which the binary layout for the protein hints should be written.
    """
    readfile = [i for i in open(inputgff, "r").readlines() if "#" not in i]
    with open(inputgff + ".code.gff", "w") as writegff:
        writegff.write("col0 \t col1 \t col2 \t col3 \t col4 \t col5 \t col6 \t col7 \t col8 \t col9\n")
        for line in readfile:
            writegff.write(line)
        writegff.close()
    read = [i.strip().split() for i in open(inputgff + ".code.gff").readlines() if i.strip().split()[2] == "mRNA"]
    mRNAplotterstart = []
    mRNAplotterend = []
    for i in range(len(read)):
        mRNAplotterstart.append(int(read[i][3]))
        mRNAplotterend.append(int(read[i][4]))
    lengthestimates = []
    for i in range(len(mRNAplotterstart)):
        lengthestimates.append(mRNAplotterend[i]-mRNAplotterstart[i])
    with open(outfile, "w") as writefile:
        writefile.write(f"These are the length estimates of the protein predicted\n")
        for i in range(len(read)):
            writefile.write(f"{mRNAplotterstart[i]}\t{mRNAplotterend[i]}\n")
        writefile.close()
    difference = []
    for i in range(len(mRNAplotterstart)):
        difference.append(mRNAplotterend[i]-mRNAplotterstart[i])
    dataframe = pd.DataFrame(difference, columns = ["mRNAlength"])
    histogram = sns.histplot(dataframe).get_figure()
    histogram.savefig("mRNAlength.png")
