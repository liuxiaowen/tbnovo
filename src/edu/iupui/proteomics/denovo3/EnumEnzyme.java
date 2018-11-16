package edu.iupui.proteomics.denovo3;

public enum EnumEnzyme {
	TRYPSIN (0, "TRYPSIN", "RK"), PEPSIN (1, "PEPSIN", "FLWYAEQ"), 
	CHYMOTRYPSIN (2, "CHYMOTRYPSIN", "FLMWY"), PROK(3, "PROK", "ACGMFSYW"), GLUC(4, "GLUC", "DE"), LYSC(5, "LYSC", "K");
	
	private int id;
	private String name;
	private String breakResidue;
	
	EnumEnzyme(int id, String name, String breakResidue) {
		this.id = id;
		this.name = name;
		this.breakResidue = breakResidue;
	}
	
	public int getId() {
		return id;
	}
	
	public String getName() {
		return name;
	}
	
	public String getBreakResidue() {
		return breakResidue;
	}
	
	public static EnumEnzyme[] getEnzymeList() {
		EnumEnzyme enzymes[] = new EnumEnzyme[6];
		enzymes[0] = EnumEnzyme.TRYPSIN;
		enzymes[1] = EnumEnzyme.PEPSIN;
		enzymes[2] = EnumEnzyme.CHYMOTRYPSIN;
		enzymes[3] = EnumEnzyme.PROK;
		enzymes[4] = EnumEnzyme.GLUC;
		enzymes[5] = EnumEnzyme.GLUC;
		
		return enzymes;
	}
}
