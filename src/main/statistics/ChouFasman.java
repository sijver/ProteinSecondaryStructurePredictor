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

    private static EnumMap<AminoAcid, EnumMap<Structure, Double>> getAminoAcidPropensities(List<Protein> proteinList){
        int[] eachAminoAcidNumber = new int[AminoAcid.values().length];
        int[] eachStructureNumber = new int[Structure.values().length];
        int[][] eachAminoAcidVsStructureNumber = new int[AminoAcid.values().length][Structure.values().length];
        int totalAminoAcidsNumber = 0;

        for(Protein protein : proteinList){
            for(int i = 0; i < protein.getProteinSequence().size(); i++){
                AminoAcid currentAminoAcid = protein.getProteinSequence().get(i);
                Structure currentStructure = protein.getProteinStructure().get(i);
                if(currentAminoAcid != null){
                    eachAminoAcidNumber[currentAminoAcid.ordinal()]++;
                    eachStructureNumber[currentStructure.ordinal()]++;
                    eachAminoAcidVsStructureNumber[currentAminoAcid.ordinal()][currentStructure.ordinal()]++;
                    totalAminoAcidsNumber++;
                }
            }
        }

        EnumMap<AminoAcid, EnumMap<Structure, Double>> aminoAcidFrequency = new EnumMap<AminoAcid, EnumMap<Structure, Double>>(AminoAcid.class);
        for(AminoAcid aminoAcid : AminoAcid.values()){
            EnumMap<Structure, Double> structurePropensites = new EnumMap<Structure, Double>(Structure.class);
            for(Structure structure : Structure.values()){
                structurePropensites.put(structure, ((double)eachAminoAcidVsStructureNumber[aminoAcid.ordinal()][structure.ordinal()] * totalAminoAcidsNumber) / ((double) eachStructureNumber[structure.ordinal()] * eachAminoAcidNumber[aminoAcid.ordinal()]));
            }
            aminoAcidFrequency.put(aminoAcid, structurePropensites);
        }


//        EnumMap<AminoAcid, EnumMap<Structure, Double>> aminoAcidFrequency = new EnumMap<AminoAcid, EnumMap<Structure, Double>>(AminoAcid.class);
//        for(AminoAcid aminoAcid : AminoAcid.values()){
//            EnumMap<Structure, Double> structurePropensites = new EnumMap<Structure, Double>(Structure.class);
//            for(Structure structure : Structure.values()){
//                structurePropensites.put(structure, 0.0);
//            }
//            aminoAcidFrequency.put(aminoAcid, structurePropensites);
//        }
//
//        for(Protein protein : proteinList){
//            for(int i = 0; i < protein.getProteinSequence().size(); i++){
//                AminoAcid currentAminoAcid = protein.getProteinSequence().get(i);
//                Structure currentStructure = protein.getProteinStructure().get(i);
//                if(currentAminoAcid != null){
//                    EnumMap<Structure, Double> structurePropensities = aminoAcidFrequency.get(currentAminoAcid);
//                    structurePropensities.put(currentStructure, structurePropensities.get(currentStructure) + 1);
//                    aminoAcidFrequency.put(currentAminoAcid, structurePropensities);
//                    totalAminoAcidsNumber++;
//                }
//            }
//        }
//
//        for(AminoAcid aminoAcid : AminoAcid.values()){
//            for(Structure structure : Structure.values()){
//                EnumMap<Structure, Double> structurePropensities = aminoAcidFrequency.get(aminoAcid);
//                structurePropensities.put(structure, structurePropensities.get(structure) / totalAminoAcidsNumber);
//                aminoAcidFrequency.put(aminoAcid, structurePropensities);
//            }
//        }

        return aminoAcidFrequency;
    }

    public static String getAminoAcidsPropensitiesString(List<Protein> proteinList){
        EnumMap<AminoAcid, EnumMap<Structure, Double>> aminoAcidPropensities = getAminoAcidPropensities(proteinList);

        StringBuilder propensitiesString = new StringBuilder();

        for(Structure structure : Structure.values()){
            propensitiesString.append("        " + structure);
        }
        propensitiesString.append("\n");

        for(AminoAcid aminoAcid : AminoAcid.values()){
            propensitiesString.append(aminoAcid);
            double maxPropensity = Double.MIN_VALUE;
            Structure maxStructure = null;
            for(Structure structure : Structure.values()){
                if(maxPropensity < aminoAcidPropensities.get(aminoAcid).get(structure)){
                    maxPropensity = aminoAcidPropensities.get(aminoAcid).get(structure);
                    maxStructure = structure;
                }
                propensitiesString.append(String.format("     %1$.4f", aminoAcidPropensities.get(aminoAcid).get(structure)));
            }
            propensitiesString.append("    (MAX = " + maxStructure + ")\n");
        }

        return propensitiesString.toString();
    }

}
