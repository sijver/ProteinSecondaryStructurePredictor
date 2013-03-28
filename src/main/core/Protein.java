package main.core;

import java.util.Arrays;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Protein {

    private List<AminoAcid> proteinSequence;
    private List<Structure> proteinStructure;
    private String pdbCode;
    private String pdbChainCode;

    public Protein(List<AminoAcid> proteinSequence, List<Structure> proteinStructure) {
        this.proteinSequence = proteinSequence;
        this.proteinStructure = proteinStructure;
    }

    public Protein(List<AminoAcid> proteinSequence, List<Structure> proteinStructure, String pdbCode, String pdbChainCode) {
        this.proteinSequence = proteinSequence;
        this.proteinStructure = proteinStructure;
        this.pdbCode = pdbCode;
        this.pdbChainCode = pdbChainCode;
    }

    public List<AminoAcid> getProteinSequence() {
        return proteinSequence;
    }

    public List<Structure> getProteinStructure() {
        return proteinStructure;
    }

    public String getPdbCode() {
        return pdbCode;
    }

    public void setPdbCode(String pdbCode) {
        this.pdbCode = pdbCode;
    }

    public String getPdbChainCode() {
        return pdbChainCode;
    }

    public void setPdbChainCode(String pdbChainCode) {
        this.pdbChainCode = pdbChainCode;
    }

    @Override
    public String toString() {
        StringBuilder proteinString = new StringBuilder();
        proteinString.append("Protein: " + pdbCode + " " + pdbChainCode + "\n");
        proteinString.append("Amino acids: " + Arrays.toString(proteinSequence.toArray()) + "\n");
        proteinString.append("Structure: " + Arrays.toString(proteinStructure.toArray()) + "\n");
        return proteinString.toString();
    }
}