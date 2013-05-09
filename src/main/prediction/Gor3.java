package main.prediction;

import main.core.AminoAcid;
import main.core.Protein;
import main.core.Structure;
import main.statistics.ChouFasman;

import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Gor3 {

    private ChouFasman chouFasman;

    private double[][] scores;

    public Gor3(ChouFasman chouFasman) {
        this.chouFasman = chouFasman;
    }

    public Protein makePrediction(Protein protein) {
        double helixScore;
        double betaScore;
        double coilScore;

        List<AminoAcid> proteinSequence = protein.getProteinSequence();

        scores = new double[proteinSequence.size()][3];

        List<Structure> predictedProteinStructure = new LinkedList<Structure>();

        for (int i = 0; i < proteinSequence.size(); i++) {
            //Self information
            helixScore = 0;
            betaScore = 0;
            coilScore = 0;
            if (proteinSequence.get(i) != null) {
                helixScore = chouFasman.getSelfInformation(proteinSequence.get(i), Structure.HELIX);
                betaScore = chouFasman.getSelfInformation(proteinSequence.get(i), Structure.BETA);
                coilScore = chouFasman.getSelfInformation(proteinSequence.get(i), Structure.COIL);


                //Pair information
                for (int j = i - ChouFasman.WINDOW_SIZE; j <= i + ChouFasman.WINDOW_SIZE; j++) {
                    if (j >= 0 && j < proteinSequence.size() && proteinSequence.get(j) != null && j != i) {
                        helixScore += chouFasman.getPairInformation(proteinSequence.get(i), proteinSequence.get(j), Structure.HELIX, j - i);
                        betaScore += chouFasman.getPairInformation(proteinSequence.get(i), proteinSequence.get(j), Structure.BETA, j - i);
                        coilScore += chouFasman.getPairInformation(proteinSequence.get(i), proteinSequence.get(j), Structure.COIL, j - i);
                    }
                }

                scores[i][0] = helixScore;
                scores[i][1] = betaScore;
                scores[i][2] = coilScore;
            } else {
                for (int j = i - ChouFasman.WINDOW_SIZE; j <= i + ChouFasman.WINDOW_SIZE; j++) {
                    if (j >= 0 && j < proteinSequence.size() && proteinSequence.get(j) != null && j != i) {
                        for (AminoAcid aminoAcid : AminoAcid.values()) {
                            helixScore += chouFasman.getPairInformation(aminoAcid, proteinSequence.get(j), Structure.HELIX, j - i) * chouFasman.getAminoAcidFrequency(aminoAcid);
                            betaScore += chouFasman.getPairInformation(aminoAcid, proteinSequence.get(j), Structure.BETA, j - i) * chouFasman.getAminoAcidFrequency(aminoAcid);
                            coilScore += chouFasman.getPairInformation(aminoAcid, proteinSequence.get(j), Structure.COIL, j - i) * chouFasman.getAminoAcidFrequency(aminoAcid);
                        }
                    }
                }
            }

            if (helixScore >= betaScore && helixScore >= coilScore) {
                predictedProteinStructure.add(Structure.HELIX);
            } else if (betaScore >= helixScore && betaScore >= coilScore) {
                predictedProteinStructure.add(Structure.BETA);
            } else {
                predictedProteinStructure.add(Structure.COIL);
            }
        }

        return new Protein(proteinSequence, predictedProteinStructure, protein.getPdbCode(), protein.getPdbChainCode());
    }

    public double[][] getScores() {
        return scores;
    }
}
