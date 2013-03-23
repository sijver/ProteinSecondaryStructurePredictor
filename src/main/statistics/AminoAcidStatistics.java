package main.statistics;

import main.core.AminoAcid;
import main.core.Protein;

import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class AminoAcidStatistics {

    private static EnumMap<AminoAcid, Double> getAminoAcidFrequency(List<Protein> proteinList){
        EnumMap<AminoAcid, Double> aminoAcidFrequency = new EnumMap<AminoAcid, Double>(AminoAcid.class);
        for(AminoAcid aminoAcid : AminoAcid.values()){
            aminoAcidFrequency.put(aminoAcid, 0.0);
        }

        int totalAminoAcidsNumber = 0;

        for(Protein protein : proteinList){
            for(AminoAcid aminoAcid : protein.getProteinSequence()){
                if(aminoAcid != null){
                    aminoAcidFrequency.put(aminoAcid, aminoAcidFrequency.get(aminoAcid) + 1);
                    totalAminoAcidsNumber++;
                }
            }
        }

        for(AminoAcid aminoAcid : AminoAcid.values()){
            aminoAcidFrequency.put(aminoAcid, aminoAcidFrequency.get(aminoAcid) / totalAminoAcidsNumber);
        }

        return aminoAcidFrequency;
    }

    public static String getAminoAcidFreqSortStats(List<Protein> proteinList){
        EnumMap<AminoAcid, Double> aminoAcidFrequency = AminoAcidStatistics.getAminoAcidFrequency(proteinList);

        StringBuilder statsString = new StringBuilder();

        for(int i = 0; i < AminoAcid.values().length; i++){
            AminoAcid maxFreqAcid = null;
            double maxFrequency = Double.MIN_VALUE;
            for(AminoAcid aminoAcid : AminoAcid.values()){
                if(aminoAcidFrequency.containsKey(aminoAcid) && maxFrequency < aminoAcidFrequency.get(aminoAcid)){
                    maxFrequency = aminoAcidFrequency.get(aminoAcid);
                    maxFreqAcid = aminoAcid;
                }
            }
            aminoAcidFrequency.remove(maxFreqAcid);
            statsString.append(maxFreqAcid + " " + maxFrequency + "\n");
        }

        return statsString.toString();
    }

    public static String getAminoAcidNameSortStats(List<Protein> proteinList){
        EnumMap<AminoAcid, Double> aminoAcidFrequency = getAminoAcidFrequency(proteinList);

        StringBuilder statsString = new StringBuilder();

        for(AminoAcid aminoAcid : AminoAcid.values()){
            statsString.append(aminoAcid + " " + aminoAcidFrequency.get(aminoAcid) + "\n");
        }
        return statsString.toString();
    }

}
