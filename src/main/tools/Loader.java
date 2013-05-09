package main.tools;

import main.core.AminoAcid;
import main.core.Protein;
import main.core.Structure;
import main.estimation.Q3;
import main.prediction.Gor3;
import main.statistics.ChouFasman;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Loader {

    private String inputFile;
    private List<Protein> proteinList;
    private List<Structure> structureList;

    public Loader(String inputFile) {
        this.inputFile = inputFile;
    }

    public List<Protein> readFile() {
        proteinList = new LinkedList<Protein>();
        BufferedReader fileReader = null;
        String line = null;

        String currentPdbCode = "";
        String currentPdbChainCode = "";
        List<AminoAcid> currentAminoAcidList = new LinkedList<AminoAcid>();
        List<Structure> currentStructureList = new LinkedList<Structure>();
        int currentPdbSequenceCode = 0;

        try {
            fileReader = new BufferedReader(new FileReader(inputFile));
            do {
                line = fileReader.readLine();
                if (line != null) {
                    String[] lineContent = line.split("[ \t]+");

                    if (!lineContent[0].equals(currentPdbCode)) {
                        if (!currentPdbCode.equals("")) {
                            proteinList.add(new Protein(currentAminoAcidList, currentStructureList, currentPdbCode, currentPdbChainCode));
                            currentAminoAcidList = new LinkedList<AminoAcid>();
                            currentStructureList = new LinkedList<Structure>();
                        }
                        currentPdbCode = lineContent[0];
                        currentPdbChainCode = lineContent[1];
                        currentPdbSequenceCode = Integer.parseInt(lineContent[2]);
                    }

                    for (int i = currentPdbSequenceCode + 1; i < Integer.parseInt(lineContent[2]); i++) {
                        currentAminoAcidList.add(null);
                        currentStructureList.add(null);
                    }

                    currentPdbSequenceCode = Integer.parseInt(lineContent[2]);

                    try {
                        currentAminoAcidList.add(AminoAcid.valueOf(lineContent[3].toUpperCase()));
                    } catch (IllegalArgumentException e) {
                        currentAminoAcidList.add(null);
                    }

                    if (lineContent[4].equalsIgnoreCase("Beta")) {
                        currentStructureList.add(Structure.BETA);
                    } else if (lineContent[4].equalsIgnoreCase("Helix")) {
                        currentStructureList.add(Structure.HELIX);
                    } else if (lineContent[4].equalsIgnoreCase("Coil") || (lineContent[4].equalsIgnoreCase("Other"))) {
                        currentStructureList.add(Structure.COIL);
                    } else {
                        currentStructureList.add(null);
                    }

                }
            } while (line != null);
        } catch (FileNotFoundException e) {
            System.out.println("File " + inputFile + " is not found");
            System.exit(1);
        } catch (IOException e) {
            System.out.println("An IO error has occured: " + e.getMessage());
            System.exit(1);
        } finally {
            proteinList.add(new Protein(currentAminoAcidList, currentStructureList, currentPdbCode, currentPdbChainCode));
            if (fileReader != null) {
                try {
                    fileReader.close();
                } catch (IOException e) {
                    // do nothing
                }
            }
        }

        return proteinList;
    }

    public List<Protein> readFastaForMSA(boolean isStructure) {
        proteinList = new LinkedList<Protein>();
        structureList = new LinkedList<Structure>();

        BufferedReader fileReader = null;
        String line = null;
        String id = "";
        StringBuffer contentBuffer = new StringBuffer();

        try {

            fileReader = new BufferedReader(new FileReader(inputFile));
            do {

                line = fileReader.readLine();
                if (line != null) {

                    line = line.trim();
                    if (line.isEmpty()) {
                        continue;
                    }
                    char firstChar = line.charAt(0);
                    if (firstChar == '>') {

                        // save the previous sequence read
                        addToProteinList(id, contentBuffer, isStructure);

                        // now can get the new id > ..
                        id = line.substring(1).trim();

                        // start a new content buffer
                        contentBuffer = new StringBuffer();

                    } else if (firstChar == ';') {

                        // comment line, skip it

                    } else {

                        // carry on reading sequence content
                        contentBuffer.append(line.trim());

                    }

                } else {

                    // save the final sequence content
                    addToProteinList(id, contentBuffer, isStructure);

                }

            } while (line != null);

        } catch (FileNotFoundException e) {
            System.out.println("File " + inputFile + " is not found");
            System.exit(1);
        } catch (IOException e) {
            System.out.println("An IO error has occured: " + e.getMessage());
            System.exit(1);
        } finally {
            if (fileReader != null) {
                try {
                    fileReader.close();
                } catch (IOException e) {
                    // do nothing
                }
            }
        }

        return proteinList;
    }

    private void addToProteinList(String id, StringBuffer sb, boolean isStructure) {
        if (sb.length() != 0) {
            if (!isStructure) {
                Protein s = null;
                String content = sb.toString();

                List<AminoAcid> aminoAcids = new LinkedList<AminoAcid>();
                for (int i = 0; i < content.length(); i++) {
                    aminoAcids.add(AminoAcid.getAminoAcidFromLetter(content.substring(i, i + 1)));
                }

                s = new Protein(aminoAcids, Collections.<Structure>emptyList(), id, "");
                this.proteinList.add(s);
            } else {
                String content = sb.toString();

                structureList = new LinkedList<Structure>();
                for (int i = 0; i < content.length(); i++) {
                    structureList.add(Structure.getStructureFromLetter(content.substring(i, i + 1)));
                }
            }
        }
    }

    public List<Structure> getStructureList() {
        return structureList;
    }

}
