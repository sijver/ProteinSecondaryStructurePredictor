package main.core;

/**
 * Created with IntelliJ IDEA.
 */
public enum Structure {
    BETA,
    HELIX,
    COIL;

    public static Structure getStructureFromLetter(String letter) {
        if (letter.equalsIgnoreCase("e")) return BETA;
        else if (letter.equalsIgnoreCase("h")) return HELIX;
        else if (letter.equalsIgnoreCase("c")) return COIL;
        else return null;
    }
}
