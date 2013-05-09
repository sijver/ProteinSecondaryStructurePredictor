package main.estimation;

import main.core.Structure;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Q3 {

    public static double getQ3Quality(List<Structure> proteinStructure1, List<Structure> proteinStructure2) {
        int residuesCorrectlyPredicted = 0;
        int residueTotal = 0;

        for (int i = 0; i < proteinStructure1.size(); i++) {
            if (proteinStructure1.get(i) == proteinStructure2.get(i)) {
                residuesCorrectlyPredicted++;
            }
            if (proteinStructure1.get(i) != null) {
                residueTotal++;
            }
        }

        return (double) residuesCorrectlyPredicted / residueTotal;
    }

}
