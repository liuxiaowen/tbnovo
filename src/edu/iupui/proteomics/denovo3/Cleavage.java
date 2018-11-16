package edu.iupui.proteomics.denovo3;

import java.util.ArrayList;

public class Cleavage {
	private double mass;
	private ArrayList<EnumEnzyme> enzymeList;
	
	public Cleavage (double mass) {
		this.mass = mass;
		enzymeList = new ArrayList<EnumEnzyme>();
	}
	
	public void addEnzyme(EnumEnzyme enzyme) {
		enzymeList.add(enzyme);
	}
	
	public double getMass() {
		return mass;
	}
	
	public ArrayList<EnumEnzyme> getEnzymeList() {
		return enzymeList;
	}

}
