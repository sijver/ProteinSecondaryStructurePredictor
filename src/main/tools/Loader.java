package main.tools;

import main.core.AminoAcid;
import main.core.Protein;
import main.core.Structure;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Loader {

    private String inputFile;
    private List<Protein> proteinList;

    public Loader(String inputFile) {
        this.inputFile = inputFile;
    }

    public List<Protein> readFile(){
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

                    if(!lineContent[0].equals(currentPdbCode)){
                        if(!currentPdbCode.equals("")){
                            proteinList.add(new Protein(currentAminoAcidList, currentStructureList, currentPdbCode, currentPdbChainCode));
                            currentAminoAcidList = new LinkedList<AminoAcid>();
                            currentStructureList = new LinkedList<Structure>();
                        }
                        currentPdbCode = lineContent[0];
                        currentPdbChainCode = lineContent[1];
                        currentPdbSequenceCode = Integer.parseInt(lineContent[2]);
                    }

                    for(int i = currentPdbSequenceCode + 1; i < Integer.parseInt(lineContent[2]); i++){
                        currentAminoAcidList.add(null);
                        currentStructureList.add(null);
                    }

                    currentPdbSequenceCode = Integer.parseInt(lineContent[2]);

                    try{
                        currentAminoAcidList.add(AminoAcid.valueOf(lineContent[3].toUpperCase()));
                    } catch (IllegalArgumentException e){
                        currentAminoAcidList.add(null);
                    }

                    if(lineContent[4].equalsIgnoreCase("Beta")){
                        currentStructureList.add(Structure.BETA);
                    } else if(lineContent[4].equalsIgnoreCase("Helix")){
                        currentStructureList.add(Structure.HELIX);
                    } else if(lineContent[4].equalsIgnoreCase("Coil") || (lineContent[4].equalsIgnoreCase("Other"))){
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
            proteinList.add(new Protein(currentAminoAcidList, currentStructureList, currentPdbCode,currentPdbChainCode));
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

}
