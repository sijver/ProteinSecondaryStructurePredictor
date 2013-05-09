package main;

import main.core.Protein;
import main.core.Structure;
import main.estimation.MatthewsCorrelationCoefficient;
import main.estimation.Q3;
import main.prediction.Gor3;
import main.statistics.ChouFasman;
import main.tools.Loader;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class StartMSA {

    public static final String[] MSA_FILE_NAMES = new String[]{"1.fasta", "2.fasta", "3.fasta", "4.fasta", "5.fasta", "6.fasta"};
    public static final String[] PROTEIN_STRUCTURE_FILE_NAMES = new String[]{"1.str", "2.str", "3.str", "4.str", "5.str", "6.str"};

    public static void main(String[] args) {
        Loader loader = new Loader("dssp_info.csv");
        List<Protein> proteinListDssp = loader.readFile();
        loader = new Loader("stride_info.csv");
        List<Protein> proteinListStride = loader.readFile();

        ChouFasman chouFasmanDssp = new ChouFasman(proteinListDssp);
        Gor3 gor3Dssp = new Gor3(chouFasmanDssp);

        ChouFasman chouFasmanStride = new ChouFasman(proteinListStride);
        Gor3 gor3Stride = new Gor3(chouFasmanStride);

        for (int f = 0; f < MSA_FILE_NAMES.length - 1; f++) {
            loader = new Loader(MSA_FILE_NAMES[f]);
            List<Protein> msaSequences = loader.readFastaForMSA(false);

            Protein predictingProtein = msaSequences.remove(0);

            loader = new Loader(PROTEIN_STRUCTURE_FILE_NAMES[f]);
            loader.readFastaForMSA(true);
            List<Structure> predictingProteinStructure = loader.getStructureList();

            List<Structure> predictingProteinStructureWithGaps = new LinkedList<Structure>();

            int counter = 0;
            for (int i = 0; i < predictingProtein.getProteinSequence().size(); i++) {
                if (predictingProtein.getProteinSequence().get(i) == null) {
                    predictingProteinStructureWithGaps.add(null);
                } else {
                    predictingProteinStructureWithGaps.add(predictingProteinStructure.get(counter));
                    counter++;
                }
            }

            double[][] scoresDssp = new double[predictingProtein.getProteinSequence().size()][3];

            double[][] scoresStride = new double[predictingProtein.getProteinSequence().size()][3];

            for (Protein msaProtein : msaSequences) {
                Protein predictedProteinDssp = gor3Dssp.makePrediction(msaProtein);
                double[][] newScores = gor3Dssp.getScores();

                for (int i = 0; i < scoresDssp.length; i++) {
                    for (int j = 0; j < 3; j++) {
                        scoresDssp[i][j] += newScores[i][j];
                    }
                }

                Protein predictedProteinStride = gor3Stride.makePrediction(msaProtein);
                newScores = gor3Stride.getScores();

                for (int i = 0; i < scoresStride.length; i++) {
                    for (int j = 0; j < 3; j++) {
                        scoresStride[i][j] += newScores[i][j];
                    }
                }
            }

            List<Structure> predictedStructureDssp = new LinkedList<Structure>();
            for (int i = 0; i < scoresDssp.length; i++) {
                if (scoresDssp[i][0] >= scoresDssp[i][1] && scoresDssp[i][0] >= scoresDssp[i][2]) {
                    predictedStructureDssp.add(Structure.HELIX);
                } else if (scoresDssp[i][1] >= scoresDssp[i][0] && scoresDssp[i][1] >= scoresDssp[i][2]) {
                    predictedStructureDssp.add(Structure.BETA);
                } else {
                    predictedStructureDssp.add(Structure.COIL);
                }
            }

            List<Structure> predictedStructureStride = new LinkedList<Structure>();
            for (int i = 0; i < scoresDssp.length; i++) {
                if (scoresStride[i][0] >= scoresStride[i][1] && scoresStride[i][0] >= scoresStride[i][2]) {
                    predictedStructureStride.add(Structure.HELIX);
                } else if (scoresStride[i][1] >= scoresStride[i][0] && scoresStride[i][1] >= scoresStride[i][2]) {
                    predictedStructureStride.add(Structure.BETA);
                } else {
                    predictedStructureStride.add(Structure.COIL);
                }
            }
            System.out.println(predictingProtein.getPdbCode() + ":");

            Protein predictedProteinDsspNoMSA = gor3Dssp.makePrediction(predictingProtein);
            Protein predictedProteinStrideNoMSA = gor3Stride.makePrediction(predictingProtein);

            System.out.println("   DSSP:");
            System.out.println("      Real structure:                    " + structureToString(predictingProteinStructure));
            System.out.println("      Predicted structure (with MSA):    " + predictStructureToString(predictingProtein, predictedStructureDssp));
            System.out.println("      Predicted structure (without MSA): " + predictStructureToString(predictingProtein, predictedProteinDsspNoMSA.getProteinStructure()));

            System.out.println(String.format("      Q3 (with/without MSA):    %1$.2f%% / %2$.2f%%", Q3.getQ3Quality(predictingProteinStructureWithGaps, predictedStructureDssp) * 100, Q3.getQ3Quality(predictingProteinStructureWithGaps, gor3Dssp.makePrediction(predictingProtein).getProteinStructure()) * 100));

            String mccHelixDsspMSA = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedStructureDssp, Structure.HELIX));
            String mccBetaDsspMSA = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedStructureDssp, Structure.BETA));
            String mccCoilDsspMSA = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedStructureDssp, Structure.COIL));

            String mccHelixDssp = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedProteinDsspNoMSA.getProteinStructure(), Structure.HELIX));
            String mccBetaDssp = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedProteinDsspNoMSA.getProteinStructure(), Structure.BETA));
            String mccCoilDssp = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedProteinDsspNoMSA.getProteinStructure(), Structure.COIL));

            System.out.println(String.format("      MCC with MSA (H/E/C):     %1$s  / %2$s  / %3$s", mccHelixDsspMSA, mccBetaDsspMSA, mccCoilDsspMSA));
            System.out.println(String.format("      MCC without MSA (H/E/C):  %1$s  / %2$s  / %3$s", mccHelixDssp, mccBetaDssp, mccCoilDssp));

            System.out.println("   STRIDE:");
            System.out.println("      Real structure:                    " + structureToString(predictingProteinStructure));
            System.out.println("      Predicted structure (with MSA):    " + predictStructureToString(predictingProtein, predictedStructureStride));
            System.out.println("      Predicted structure (without MSA): " + predictStructureToString(predictingProtein, predictedProteinStrideNoMSA.getProteinStructure()));

            System.out.println(String.format("      Q3 (with/without MSA):    %1$.2f%% / %2$.2f%%", Q3.getQ3Quality(predictingProteinStructureWithGaps, predictedStructureStride) * 100, Q3.getQ3Quality(predictingProteinStructureWithGaps, gor3Stride.makePrediction(predictingProtein).getProteinStructure()) * 100));

            String mccHelixStrideMSA = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedStructureStride, Structure.HELIX));
            String mccBetaStrideMSA = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedStructureStride, Structure.BETA));
            String mccCoilStrideMSA = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedStructureStride, Structure.COIL));

            String mccHelixStride = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedProteinStrideNoMSA.getProteinStructure(), Structure.HELIX));
            String mccBetaStride = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedProteinStrideNoMSA.getProteinStructure(), Structure.BETA));
            String mccCoilStride = String.format("%1$.3f", MatthewsCorrelationCoefficient.getMatthewsCorrelationCoefficient(predictingProteinStructureWithGaps, predictedProteinStrideNoMSA.getProteinStructure(), Structure.COIL));

            System.out.println(String.format("      MCC with MSA (H/E/C):     %1$s  / %2$s  / %3$s", mccHelixStrideMSA, mccBetaStrideMSA, mccCoilStrideMSA));
            System.out.println(String.format("      MCC without MSA (H/E/C):  %1$s  / %2$s  / %3$s", mccHelixStride, mccBetaStride, mccCoilStride));

            System.out.println();
        }

    }

    public static String structureToString(List<Structure> structureList) {
        StringBuilder string = new StringBuilder();

        for (Structure structure : structureList) {
            string.append(structure.toString().substring(0, 1));
        }
        return string.toString();
    }

    public static String predictStructureToString(Protein realProtein, List<Structure> predictingStructureList) {
        StringBuilder stringPredict = new StringBuilder();

        for (int i = 0; i < realProtein.getProteinSequence().size(); i++) {
            if (realProtein.getProteinSequence().get(i) != null) {
                stringPredict.append(predictingStructureList.get(i).toString().substring(0, 1));
            }
        }
        return stringPredict.toString();
    }

    public static List<Structure> getStructureWithoutGaps(Protein predictingProteinWithGaps, List<Structure> structureList) {
        List<Structure> structure = new ArrayList<Structure>();
        for (int i = 0; i < predictingProteinWithGaps.getProteinSequence().size(); i++) {
            if (predictingProteinWithGaps.getProteinSequence().get(i) != null) {
                structure.add(structureList.get(i));
            }
        }
        return structure;
    }

}
