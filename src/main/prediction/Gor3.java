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

    ChouFasman chouFasman;

    public Gor3(ChouFasman chouFasman) {
        this.chouFasman = chouFasman;
    }

    public Protein makePrediction(Protein protein){
        double helixScore;
        double betaScore;
        double coilScore;
        int windowSize = 8;

        List<AminoAcid> proteinSequence = protein.getProteinSequence();

        List<Structure> predictedProteinStructure = new LinkedList<Structure>();

        for(int i = 0; i < proteinSequence.size(); i++){
            //Self information
            helixScore = chouFasman.getSelfInformation(proteinSequence.get(i), Structure.HELIX);
            betaScore = chouFasman.getSelfInformation(proteinSequence.get(i), Structure.BETA);
            coilScore = chouFasman.getSelfInformation(proteinSequence.get(i), Structure.COIL);

            //Pair information
            for(int j = i - windowSize; j <= i + windowSize; j++){
                if(j >= 0 && j < proteinSequence.size() && proteinSequence.get(j) != null && j != i){
                    helixScore += chouFasman.getPairInformation(proteinSequence.get(i), proteinSequence.get(j), Structure.HELIX);
                    betaScore += chouFasman.getPairInformation(proteinSequence.get(i), proteinSequence.get(j), Structure.BETA);
                    coilScore += chouFasman.getPairInformation(proteinSequence.get(i), proteinSequence.get(j), Structure.COIL);
                }
            }
            if(helixScore >= betaScore && helixScore >= coilScore){
                predictedProteinStructure.add(Structure.HELIX);
            } else if(betaScore >= helixScore && betaScore >= coilScore){
                predictedProteinStructure.add(Structure.BETA);
            } else {
                predictedProteinStructure.add(Structure.COIL);
            }
        }

        return new Protein(proteinSequence, predictedProteinStructure, protein.getPdbCode(), protein.getPdbChainCode());
    }

}
