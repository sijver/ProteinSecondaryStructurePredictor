package main.core;

/**
 * Created with IntelliJ IDEA.
 */
public enum AminoAcid {
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLN,
    GLU,
    GLY,
    HIS,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL;

    public static AminoAcid getAminoAcidFromLetter(String letter) {
        if (letter.equalsIgnoreCase("a")) return ALA;
        else if (letter.equalsIgnoreCase("r")) return ARG;
        else if (letter.equalsIgnoreCase("n")) return ASN;
        else if (letter.equalsIgnoreCase("d")) return ASP;
        else if (letter.equalsIgnoreCase("c")) return CYS;
        else if (letter.equalsIgnoreCase("q")) return GLN;
        else if (letter.equalsIgnoreCase("e")) return GLU;
        else if (letter.equalsIgnoreCase("g")) return GLY;
        else if (letter.equalsIgnoreCase("h")) return HIS;
        else if (letter.equalsIgnoreCase("i")) return ILE;
        else if (letter.equalsIgnoreCase("l")) return LEU;
        else if (letter.equalsIgnoreCase("k")) return LYS;
        else if (letter.equalsIgnoreCase("m")) return MET;
        else if (letter.equalsIgnoreCase("f")) return PHE;
        else if (letter.equalsIgnoreCase("p")) return PRO;
        else if (letter.equalsIgnoreCase("s")) return SER;
        else if (letter.equalsIgnoreCase("t")) return THR;
        else if (letter.equalsIgnoreCase("w")) return TRP;
        else if (letter.equalsIgnoreCase("y")) return TYR;
        else if (letter.equalsIgnoreCase("v")) return VAL;
        else return null;
    }
}
