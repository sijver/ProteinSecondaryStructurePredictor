package main.statistics;

import main.core.AminoAcid;
import main.core.Protein;
import main.core.Structure;

import java.util.EnumMap;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class ChouFasman {

    private EnumMap<AminoAcid, EnumMap<Structure, Double>> aminoAcidPropensities = null;
    private double[][] selfInformation;
    private double[][][] pairInformation;

    public ChouFasman(List<Protein> proteinList) {
        computeAminoAcidPropensities(proteinList);
    }

    private void computeAminoAcidPropensities(List<Protein> proteinList) {
        int[] eachAminoAcidNumber = new int[AminoAcid.values().length];
        int[] eachStructureNumber = new int[Structure.values().length];
        int[][] eachAminoAcidVsStructureNumber = new int[AminoAcid.values().length][Structure.values().length];
        int totalAminoAcidsNumber = 0;

        for (Protein protein : proteinList) {
            for (int i = 0; i < protein.getProteinSequence().size(); i++) {
                AminoAcid currentAminoAcid = protein.getProteinSequence().get(i);
                Structure currentStructure = protein.getProteinStructure().get(i);
                if (currentAminoAcid != null) {
                    eachAminoAcidNumber[currentAminoAcid.ordinal()]++;
                    eachStructureNumber[currentStructure.ordinal()]++;
                    eachAminoAcidVsStructureNumber[currentAminoAcid.ordinal()][currentStructure.ordinal()]++;
                    totalAminoAcidsNumber++;
                }
            }
        }

        aminoAcidPropensities = new EnumMap<AminoAcid, EnumMap<Structure, Double>>(AminoAcid.class);
        for (AminoAcid aminoAcid : AminoAcid.values()) {
            EnumMap<Structure, Double> structurePropensites = new EnumMap<Structure, Double>(Structure.class);
            for (Structure structure : Structure.values()) {
                structurePropensites.put(structure, ((double) eachAminoAcidVsStructureNumber[aminoAcid.ordinal()][structure.ordinal()] * totalAminoAcidsNumber) / ((double) eachStructureNumber[structure.ordinal()] * eachAminoAcidNumber[aminoAcid.ordinal()]));
            }
            aminoAcidPropensities.put(aminoAcid, structurePropensites);
        }

        selfInformation = new double[AminoAcid.values().length][Structure.values().length];

        for (int i = 0; i < AminoAcid.values().length; i++) {
            for (int j = 0; j < Structure.values().length; j++) {
                selfInformation[i][j] = Math.log((double) eachAminoAcidVsStructureNumber[i][j] / (eachAminoAcidNumber[i] - eachAminoAcidVsStructureNumber[i][j])) + Math.log(((double) totalAminoAcidsNumber - eachStructureNumber[j]) / eachStructureNumber[j]);
            }
        }

        pairInformation = new double[AminoAcid.values().length][AminoAcid.values().length][Structure.values().length];

        for (int i = 0; i < AminoAcid.values().length; i++) {
            for (int j = 0; j < AminoAcid.values().length; j++) {
                for (int k = 0; k < Structure.values().length; k++) {
                    pairInformation[i][j][k] = Math.log((double) eachAminoAcidVsStructureNumber[i][j] / (eachAminoAcidNumber[i] - eachAminoAcidVsStructureNumber[i][j])) + Math.log(((double) totalAminoAcidsNumber - eachStructureNumber[j]) / eachStructureNumber[j]);
                }
            }
        }

    }

    public String getAminoAcidsPropensitiesString() {
        StringBuilder propensitiesString = new StringBuilder();

        for (Structure structure : Structure.values()) {
            propensitiesString.append("        " + structure);
        }
        propensitiesString.append("\n");

        for (AminoAcid aminoAcid : AminoAcid.values()) {
            propensitiesString.append(aminoAcid);
            double maxPropensity = Double.MIN_VALUE;
            Structure maxStructure = null;
            for (Structure structure : Structure.values()) {
                if (maxPropensity < aminoAcidPropensities.get(aminoAcid).get(structure)) {
                    maxPropensity = aminoAcidPropensities.get(aminoAcid).get(structure);
                    maxStructure = structure;
                }
                propensitiesString.append(String.format("     %1$.4f", aminoAcidPropensities.get(aminoAcid).get(structure)));
            }
            propensitiesString.append("    (MAX = " + maxStructure + ")\n");
        }

        return propensitiesString.toString();
    }

    public double getSelfInformation(AminoAcid aminoAcid, Structure structure) {
        return selfInformation[aminoAcid.ordinal()][structure.ordinal()];
    }

    public double getPairInformation(AminoAcid aminoAcid1, AminoAcid aminoAcid2, Structure structure) {
        return pairInformation[aminoAcid1.ordinal()][aminoAcid2.ordinal()][structure.ordinal()];
    }

}
