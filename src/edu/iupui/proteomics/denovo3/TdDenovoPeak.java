package edu.iupui.proteomics.denovo3;

import edu.iupui.proteomics.spec.deconvsp.DeconvPeak;
import edu.iupui.proteomics.spec.peak.Peak;
import edu.iupui.proteomics.spec.prmsp.EnumPrmPeakType;

public class TdDenovoPeak implements Peak{
	private DeconvPeak basePeak; /* base vertex */
	private double monoMass;
	private double intensity;
	protected double score = 1;
	private EnumPrmPeakType baseType; 
	private double tolerance;

	public TdDenovoPeak(DeconvPeak basePeak, double monoMass,
			EnumPrmPeakType baseType, double score) throws Exception {
		this.monoMass = monoMass;
		this.intensity = basePeak.getIntensity();
		this.basePeak = basePeak;
		this.score = score;
		this.baseType = baseType;
	}

	/* gets */
	public DeconvPeak getBasePeak() {
		return basePeak;
	}

	public double getIntensity() {
		return intensity;
	}

	public double getMonoMass() {
		return monoMass;
	}

	public double getPosition() {
		return monoMass;
	}

	public double getTolerance() {
		return tolerance;
	}

	public void setPosition(double p) {
		monoMass = p;
	}

	public void setIntensity(double i) {
		intensity = i;
	}

	public void setTolerance(double tolerance) {
		this.tolerance = tolerance;
	}

	public EnumPrmPeakType getBaseType() {
		return baseType;
	}


	public double getScore() {
		return score;
	}

	public void increaseScore() {
		score = score + 1.0;
	}
}
