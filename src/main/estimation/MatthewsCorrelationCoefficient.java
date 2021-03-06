package main.estimation;

import main.core.Structure;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class MatthewsCorrelationCoefficient {

    public static double getMatthewsCorrelationCoefficient(List<Structure> proteinStructure1, List<Structure> proteinStructure2, Structure structure) {
        int truePositive = 0;
        int trueNegative = 0;
        int falsePositive = 0;
        int falseNegative = 0;

        for (int i = 0; i < proteinStructure1.size(); i++) {
            if (proteinStructure1.get(i) != null && proteinStructure2.get(i) != null) {
                if (proteinStructure1.get(i) == structure) {
                    if (proteinStructure2.get(i) == structure) {
                        truePositive++;
                    } else {
                        falseNegative++;
                    }
                } else {
                    if (proteinStructure2.get(i) == structure) {
                        falsePositive++;
                    } else if (proteinStructure1.get(i) == proteinStructure2.get(i)) {
                        trueNegative++;
                    }
                }
            }
        }

//        if (truePositive + falsePositive == 0 || truePositive + falseNegative == 0 || trueNegative + falsePositive == 0 || trueNegative + falseNegative == 0) {
//            return 0;
//        } else {
//            return ((double) truePositive * trueNegative - falsePositive * falseNegative) / Math.sqrt((truePositive + falsePositive) * (truePositive + falseNegative) * (trueNegative + falsePositive) * (trueNegative + falseNegative));
//        }
        if (truePositive == 0)
            truePositive = 1;
        if (trueNegative == 0)
            trueNegative = 1;
        if (falsePositive == 0)
            falsePositive = 1;
        if (falseNegative == 0)
            falseNegative = 1;

        return ((double) truePositive * trueNegative - falsePositive * falseNegative) / Math.sqrt((truePositive + falsePositive) * (truePositive + falseNegative) * (trueNegative + falsePositive) * (trueNegative + falseNegative));

    }

}
