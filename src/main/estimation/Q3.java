package main.estimation;

import main.core.Structure;

import java.util.List;

/**
 * Created with IntelliJ IDEA.
 */
public class Q3 {

    public double getQ3Quality(List<Structure> proteinStructure1, List<Structure> proteinStructure2){
        int residuesCorrectlyPredicted = 0;

        for(int i = 0; i < proteinStructure1.size(); i++){
            if(proteinStructure1.get(i) == proteinStructure2.get(i)){
                residuesCorrectlyPredicted++;
            }
        }

        return (double)residuesCorrectlyPredicted / proteinStructure1.size();
    }

}
